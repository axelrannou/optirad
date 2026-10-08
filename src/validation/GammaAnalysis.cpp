#include "validation/GammaAnalysis.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace optirad {

namespace {

struct SearchOffset {
    double r2;          // squared distance [mm^2]
    Vec3   voxelDelta;  // offset expressed in evaluated-grid voxel coordinates
};

// Returns false if the point is outside the evaluated grid (no extrapolation).
bool sampleTrilinear(const double* data, const std::array<size_t, 3>& dims,
                     double fi, double fj, double fk, double& out) {
    const double nx = static_cast<double>(dims[0]);
    const double ny = static_cast<double>(dims[1]);
    const double nz = static_cast<double>(dims[2]);
    if (fi < 0.0 || fj < 0.0 || fk < 0.0 || fi > nx - 1.0 || fj > ny - 1.0 || fk > nz - 1.0)
        return false;

    const size_t i0 = static_cast<size_t>(fi), j0 = static_cast<size_t>(fj), k0 = static_cast<size_t>(fk);
    const size_t i1 = std::min(i0 + 1, dims[0] - 1);
    const size_t j1 = std::min(j0 + 1, dims[1] - 1);
    const size_t k1 = std::min(k0 + 1, dims[2] - 1);
    const double di = fi - static_cast<double>(i0);
    const double dj = fj - static_cast<double>(j0);
    const double dk = fk - static_cast<double>(k0);

    auto at = [&](size_t i, size_t j, size_t k) { return data[i + j * dims[0] + k * dims[0] * dims[1]]; };
    const double c00 = at(i0, j0, k0) * (1 - di) + at(i1, j0, k0) * di;
    const double c10 = at(i0, j1, k0) * (1 - di) + at(i1, j1, k0) * di;
    const double c01 = at(i0, j0, k1) * (1 - di) + at(i1, j0, k1) * di;
    const double c11 = at(i0, j1, k1) * (1 - di) + at(i1, j1, k1) * di;
    out = (c00 * (1 - dj) + c10 * dj) * (1 - dk) + (c01 * (1 - dj) + c11 * dj) * dk;
    return true;
}

} // namespace

GammaResult GammaAnalysis::compute(const DoseMatrix& reference, const DoseMatrix& evaluated,
                                   const GammaOptions& opt) {
    if (opt.doseDifferencePercent <= 0.0 || opt.distanceToAgreementMm <= 0.0 ||
        opt.searchStepFraction <= 0.0 || opt.searchRadiusFactor <= 0.0)
        throw std::invalid_argument("GammaAnalysis: criteria and search parameters must be positive");
    if (reference.size() == 0 || evaluated.size() == 0)
        throw std::invalid_argument("GammaAnalysis: empty dose matrix");

    const Grid& refGrid  = reference.getGrid();
    const Grid& evalGrid = evaluated.getGrid();
    const auto refDims   = refGrid.getDimensions();
    const auto evalDims  = evalGrid.getDimensions();

    const double dta = opt.distanceToAgreementMm;
    const double maxGamma = opt.searchRadiusFactor;
    const double step = opt.searchStepFraction * dta;
    const double radius = maxGamma * dta;
    const int n = static_cast<int>(std::ceil(radius / step));

    // Offsets are precomputed serially: Grid caches its matrices lazily and is not thread-safe.
    const Vec3 evalOriginVoxel = evalGrid.patientToVoxel(evalGrid.getOrigin());
    std::vector<SearchOffset> offsets;
    for (int a = -n; a <= n; ++a)
        for (int b = -n; b <= n; ++b)
            for (int c = -n; c <= n; ++c) {
                const double ox = a * step, oy = b * step, oz = c * step;
                const double r2 = ox * ox + oy * oy + oz * oz;
                if (r2 > radius * radius) continue;
                const Vec3 p = vecAdd(evalGrid.getOrigin(), Vec3{ox, oy, oz});
                offsets.push_back({r2, vecSub(evalGrid.patientToVoxel(p), evalOriginVoxel)});
            }
    std::sort(offsets.begin(), offsets.end(),
              [](const SearchOffset& x, const SearchOffset& y) { return x.r2 < y.r2; });

    const double refMax = reference.getMax();
    const double cutoff = refMax * opt.lowDoseCutoffPercent / 100.0;
    const double ddFrac = opt.doseDifferencePercent / 100.0;
    const double ddGlobal = ddFrac * refMax;

    const size_t nx = refDims[0], ny = refDims[1], nz = refDims[2];
    GammaResult res;
    res.gamma.assign(reference.size(), std::numeric_limits<double>::quiet_NaN());

    // Voxel centres of the reference grid mapped into evaluated-grid voxel space.
    // Each is computed serially here for the same thread-safety reason.
    std::vector<Vec3> refToEval(reference.size());
    for (size_t k = 0; k < nz; ++k)
        for (size_t j = 0; j < ny; ++j)
            for (size_t i = 0; i < nx; ++i) {
                const size_t idx = i + j * nx + k * nx * ny;
                if (reference.data()[idx] < cutoff || reference.data()[idx] <= 0.0) continue;
                refToEval[idx] = evalGrid.patientToVoxel(
                    refGrid.voxelToPatient(Vec3{double(i), double(j), double(k)}));
            }

    const double* refData  = reference.data();
    const double* evalData = evaluated.data();

#pragma omp parallel for schedule(dynamic, 64)
    for (long long idx = 0; idx < static_cast<long long>(reference.size()); ++idx) {
        const double dRef = refData[idx];
        if (dRef < cutoff || dRef <= 0.0) continue;

        const double dd = opt.globalNormalization ? ddGlobal : ddFrac * dRef;
        const Vec3& c = refToEval[idx];
        double best2 = maxGamma * maxGamma;

        for (const auto& off : offsets) {
            const double spatial2 = off.r2 / (dta * dta);
            if (spatial2 >= best2) break;  // offsets are sorted: nothing closer remains
            double dEval;
            if (!sampleTrilinear(evalData, evalDims, c[0] + off.voxelDelta[0],
                                 c[1] + off.voxelDelta[1], c[2] + off.voxelDelta[2], dEval))
                continue;
            const double diff = (dEval - dRef) / dd;
            best2 = std::min(best2, spatial2 + diff * diff);
        }
        res.gamma[idx] = std::sqrt(best2);
    }

    std::vector<double> valid;
    for (double g : res.gamma)
        if (!std::isnan(g)) valid.push_back(g);

    res.numEvaluated = valid.size();
    if (!valid.empty()) {
        double sum = 0.0;
        for (double g : valid) {
            sum += g;
            if (g <= 1.0) ++res.numPassed;
        }
        res.meanGamma = sum / valid.size();
        res.passRate  = 100.0 * res.numPassed / valid.size();
        std::sort(valid.begin(), valid.end());
        res.maxGamma = valid.back();
        res.gamma95  = valid[std::min(valid.size() - 1, static_cast<size_t>(0.95 * valid.size()))];
    }
    return res;
}

} // namespace optirad
