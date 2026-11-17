#include "NormalTissueObjective.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>

namespace optirad {

namespace {

constexpr double kFar = 1e20;   // squared distance of "no target voxel yet"

// Exact 1D squared distance transform (Felzenszwalb & Huttenlocher) of samples `spacing` mm apart.
void distanceTransform1D(std::vector<double>& f, double spacing) {
    const size_t n = f.size();
    if (n < 2) return;
    std::vector<size_t> v(n);
    std::vector<double> z(n + 1), out(n);
    auto x = [spacing](size_t i) { return spacing * static_cast<double>(i); };
    auto intersection = [&](size_t q, size_t p) {
        return ((f[q] + x(q) * x(q)) - (f[p] + x(p) * x(p))) / (2.0 * (x(q) - x(p)));
    };
    size_t k = 0;
    z[0] = -std::numeric_limits<double>::infinity();
    z[1] = std::numeric_limits<double>::infinity();
    for (size_t q = 1; q < n; ++q) {
        double s = intersection(q, v[k]);
        while (s <= z[k]) {
            --k;
            s = intersection(q, v[k]);
        }
        ++k;
        v[k] = q;
        z[k] = s;
        z[k + 1] = std::numeric_limits<double>::infinity();
    }
    k = 0;
    for (size_t q = 0; q < n; ++q) {
        while (z[k + 1] < x(q)) ++k;
        const double d = x(q) - x(v[k]);
        out[q] = std::min(kFar, d * d + f[v[k]]);
    }
    f = std::move(out);
}

// Euclidean distance (mm) of every voxel to the nearest voxel of `inside`.
std::vector<double> distanceToMask(const std::vector<bool>& inside, const Grid& grid) {
    const auto dims = grid.getDimensions();
    const auto spacing = grid.getSpacing();
    const long long nx = static_cast<long long>(dims[0]), ny = static_cast<long long>(dims[1]),
                    nz = static_cast<long long>(dims[2]);
    std::vector<double> d(inside.size());
    for (size_t i = 0; i < inside.size(); ++i) d[i] = inside[i] ? 0.0 : kFar;
    // One pass per axis: `stride` between samples of a line, lines enumerated over the two other axes.
    auto pass = [&](long long length, double step, long long stride, long long countA, long long strideA,
                    long long countB, long long strideB) {
        #pragma omp parallel for collapse(2)
        for (long long b = 0; b < countB; ++b) {
            for (long long a = 0; a < countA; ++a) {
                std::vector<double> line(static_cast<size_t>(length));
                const long long base = a * strideA + b * strideB;
                for (long long i = 0; i < length; ++i) line[i] = d[base + i * stride];
                distanceTransform1D(line, step);
                for (long long i = 0; i < length; ++i) d[base + i * stride] = line[i];
            }
        }
    };
    pass(nx, spacing[0], 1, ny, nx, nz, nx * ny);
    pass(ny, spacing[1], nx, nx, 1, nz, nx * ny);
    pass(nz, spacing[2], nx * ny, nx, 1, ny, nx);
    for (double& value : d) value = std::sqrt(value);
    return d;
}

} // namespace

NormalTissueObjective::NormalTissueObjective(std::vector<size_t> voxels, std::vector<double> limits,
                                             std::string summary)
    : m_limits(std::move(limits)), m_summary(std::move(summary)) {
    setVoxelIndices(voxels);
}

double NormalTissueObjective::limit(const NtoSettings& settings, double prescriptionDose, double distanceMm) {
    const double start = 0.01 * settings.startDosePct * prescriptionDose;
    const double end = 0.01 * settings.endDosePct * prescriptionDose;
    const double beyond = distanceMm - settings.distanceFromTargetMm;
    if (beyond <= 0.0) return start;
    const double decay = std::exp(-settings.fallOffPerMm * beyond);
    return start * decay + end * (1.0 - decay);
}

std::unique_ptr<NormalTissueObjective> NormalTissueObjective::build(
    const NtoSettings& settings, double prescriptionDose, const StructureSet* structures, const Grid& ctGrid,
    const Grid& doseGrid, const IDoseOperator* dose) {
    if (!settings.enabled || settings.priority <= 0.0 || prescriptionDose <= 0.0) return nullptr;
    if (!structures) return nullptr;
    const size_t voxels = doseGrid.getNumVoxels();
    std::vector<bool> target(voxels, false), body(voxels, false);
    bool hasTarget = false, hasBody = false;
    for (size_t i = 0; i < structures->getCount(); ++i) {
        const auto* structure = structures->getStructure(i);
        if (!structure || structure->getVoxelIndices().empty()) continue;
        if (!structure->isTarget() && !structure->isExternal()) continue;
        auto& mask = structure->isTarget() ? target : body;
        for (size_t voxel : Grid::mapVoxelIndices(ctGrid, doseGrid, structure->getVoxelIndices()))
            if (voxel < voxels) mask[voxel] = true;
        (structure->isTarget() ? hasTarget : hasBody) = true;
    }
    if (!hasTarget) {
        Logger::warn("NTO: no target structure; normal tissue objective disabled");
        return nullptr;
    }
    if (!hasBody) Logger::warn("NTO: no BODY/EXTERNAL structure; applying the NTO to the whole dose grid");

    // Voxel centres: the target border lies half a voxel from the centre of the outermost target voxel.
    const auto spacing = doseGrid.getSpacing();
    const double border = 0.5 * std::min({spacing[0], spacing[1], spacing[2]});
    const auto distance = distanceToMask(target, doseGrid);
    std::vector<size_t> indices;
    std::vector<double> limits;
    for (size_t voxel = 0; voxel < voxels; ++voxel) {
        if (target[voxel] || (hasBody && !body[voxel])) continue;
        if (dose && (voxel >= dose->getNumVoxels() || !dose->isRowComputed(voxel))) continue;
        indices.push_back(voxel);
        limits.push_back(limit(settings, prescriptionDose, std::max(0.0, distance[voxel] - border)));
    }
    std::ostringstream summary;
    summary << "start " << settings.startDosePct << "% up to " << settings.distanceFromTargetMm
            << " mm, end " << settings.endDosePct << "%, fall-off " << settings.fallOffPerMm << "/mm";
    std::ostringstream curve;
    curve << std::fixed << std::setprecision(1) << "NTO limit vs distance from target border:";
    for (double mm : {0.0, 5.0, 10.0, 20.0, 30.0, 50.0})
        curve << ' ' << mm << " mm=" << limit(settings, prescriptionDose, mm) << " Gy";
    Logger::info(curve.str());
    auto objective = std::make_unique<NormalTissueObjective>(std::move(indices), std::move(limits), summary.str());
    objective->setWeight(settings.priority);
    return objective;
}

std::string NormalTissueObjective::describe() const {
    std::ostringstream out;
    out << "NTO (" << m_summary << ") priority=" << m_weight << " voxels=" << getActiveIndices().size();
    return out.str();
}

double NormalTissueObjective::compute(const std::vector<double>& dose) const {
    const auto& indices = getActiveIndices();
    if (indices.empty()) return 0.0;
    double sum = 0.0;
    #pragma omp parallel for reduction(+:sum)
    for (long long i = 0; i < static_cast<long long>(indices.size()); ++i) {
        if (indices[i] >= dose.size()) continue;
        const double excess = dose[indices[i]] - m_limits[i];
        if (excess > 0.0) sum += excess * excess;
    }
    return m_weight * sum / static_cast<double>(indices.size());
}

std::vector<double> NormalTissueObjective::gradient(const std::vector<double>& dose) const {
    std::vector<double> grad(dose.size(), 0.0);
    const auto& indices = getActiveIndices();
    if (indices.empty()) return grad;
    const double scale = 2.0 * m_weight / static_cast<double>(indices.size());
    #pragma omp parallel for
    for (long long i = 0; i < static_cast<long long>(indices.size()); ++i) {
        if (indices[i] >= dose.size()) continue;
        const double excess = dose[indices[i]] - m_limits[i];
        if (excess > 0.0) grad[indices[i]] = scale * excess;
    }
    return grad;
}

} // namespace optirad
