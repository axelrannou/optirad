#include "validation/WaterPhantom.hpp"

#include "core/workflow/DoseCalculationPipeline.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "core/Beam.hpp"
#include "core/Ray.hpp"
#include "dose/DoseEngineFactory.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>

namespace optirad {

WaterPhantomResult WaterPhantom::compute(const WaterPhantomConfig& cfg) {
    const double res = cfg.resolutionMm;
    const size_t nx = size_t(std::round(2.0 * cfg.lateralHalfSizeMm / res));
    const size_t ny = size_t(std::round(cfg.depthMm / res));
    const size_t nz = nx;
    const Vec3 origin = {-cfg.lateralHalfSizeMm, -cfg.isocenterDepthMm, -cfg.lateralHalfSizeMm};

    // Grid index order is [y, x, z] (first index runs along patient y, as in matRad).
    Grid grid;
    grid.setDimensions(ny, nx, nz);
    grid.setSpacing(res, res, res);
    grid.setOrigin(origin);

    auto ct = std::make_unique<Volume<int16_t>>();
    ct->setGrid(grid);
    ct->allocate();
    for (size_t i = 0; i < ct->size(); ++i) ct->data()[i] = 0;  // 0 HU = water

    // Target box centred on the isocenter, cross-section equal to the field size.
    const double half = cfg.fieldSizeMm / 2.0;
    std::vector<size_t> targetVoxels;
    for (size_t k = 0; k < nz; ++k)
        for (size_t j = 0; j < nx; ++j)
            for (size_t i = 0; i < ny; ++i) {
                const double x = origin[0] + j * res, y = origin[1] + i * res, z = origin[2] + k * res;
                if (std::abs(x) <= half && std::abs(z) <= half && std::abs(y) <= 2.0 * res)
                    targetVoxels.push_back(i + j * ny + k * ny * nx);
            }
    if (targetVoxels.empty()) throw std::runtime_error("WaterPhantom: empty target");

    auto ptv = std::make_unique<Structure>();
    ptv->setName("PTV");
    ptv->setType("PTV");
    ptv->setVoxelIndices(targetVoxels);
    ptv->setPreRasterized(true);

    // The engine only scores voxels inside structures, so a body covering the whole phantom is required.
    std::vector<size_t> bodyVoxels(nx * ny * nz);
    for (size_t i = 0; i < bodyVoxels.size(); ++i) bodyVoxels[i] = i;
    auto body = std::make_unique<Structure>();
    body->setName("BODY");
    body->setType("EXTERNAL");
    body->setVoxelIndices(bodyVoxels);
    body->setPreRasterized(true);

    auto structures = std::make_unique<StructureSet>();
    structures->addStructure(std::move(ptv));
    structures->addStructure(std::move(body));

    auto patient = std::make_shared<PatientData>();
    patient->setCTVolume(std::move(ct));
    patient->setStructureSet(std::move(structures));
    patient->convertHUtoED();

    PlanConfig pc;
    pc.machineName = cfg.machineName;
    pc.bixelWidth = cfg.bixelWidthMm;
    pc.gantryAngles = {0.0};
    auto built = PlanBuilder::build(pc, patient);
    if (!built.stf) throw std::runtime_error("WaterPhantom: STF generation failed");

    DoseCalcPipelineOptions opts;
    opts.resolution = {res, res, res};
    opts.useCache = false;
    opts.relativeThreshold = cfg.relativeThreshold;
    opts.absoluteThreshold = 0.0;
    const auto t0 = std::chrono::steady_clock::now();
    auto dijRes = DoseCalculationPipeline::run(*built.plan, *built.stf, *patient, opts);
    const auto t1 = std::chrono::steady_clock::now();

    auto engine = DoseEngineFactory::create("PencilBeam");
    const size_t nBixels = dijRes.dij->getNumBixels();
    const Beam* beam = built.stf->getBeam(0);
    if (beam->getNumOfRays() != nBixels) throw std::runtime_error("WaterPhantom: expected one bixel per ray");

    WaterPhantomResult out;
    const size_t nDoses = std::max<size_t>(1, cfg.fluences.size());
    for (size_t f = 0; f < nDoses; ++f) {
        std::vector<double> weights(nBixels, 1.0);
        if (!cfg.fluences.empty()) {
            for (size_t r = 0; r < nBixels; ++r) {
                const Vec3 p = beam->getRay(r)->getRayPosBev();
                weights[r] = cfg.fluences[f](p[0], p[2]);
            }
        }
        out.doses.push_back(std::make_shared<DoseMatrix>(
            engine->calculateDose(*dijRes.dij, weights, *dijRes.doseGrid)));
    }
    out.dose = out.doses.front();
    out.dijNonZeros = dijRes.dij->getNumNonZeros();
    out.dijSeconds = std::chrono::duration<double>(t1 - t0).count();
    out.surfaceY = origin[1];
    return out;
}

namespace {

// Trilinear sample; positions outside the grid are skipped.
bool sampleAt(const DoseMatrix& dose, double x, double y, double z, double& out) {
    const Vec3 v = dose.getGrid().patientToVoxel(Vec3{x, y, z});
    const auto d = dose.getGrid().getDimensions();
    if (v[0] < 0 || v[1] < 0 || v[2] < 0 || v[0] > d[0] - 1.0 || v[1] > d[1] - 1.0 || v[2] > d[2] - 1.0)
        return false;
    out = dose.interpolateAt(v[0], v[1], v[2]);
    return true;
}

} // namespace

DoseCurve sampleDepthDose(const DoseMatrix& dose, double x, double z,
                          double yStart, double yEnd, double stepMm) {
    DoseCurve c;
    for (double y = yStart; y <= yEnd + 1e-9; y += stepMm) {
        double d;
        if (sampleAt(dose, x, y, z, d)) { c.position.push_back(y); c.dose.push_back(d); }
    }
    return c;
}

DoseCurve sampleProfileX(const DoseMatrix& dose, double y, double z,
                         double xStart, double xEnd, double stepMm) {
    DoseCurve c;
    for (double x = xStart; x <= xEnd + 1e-9; x += stepMm) {
        double d;
        if (sampleAt(dose, x, y, z, d)) { c.position.push_back(x); c.dose.push_back(d); }
    }
    return c;
}

} // namespace optirad
