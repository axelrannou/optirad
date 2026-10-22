#pragma once

#include "GridInterpolant2D.hpp"
#include "core/Beam.hpp"
#include <memory>
#include <vector>

namespace optirad {

class CudaDirectDose {
public:
    explicit CudaDirectDose(size_t numVoxels);
    ~CudaDirectDose();
    void accumulate(const Beam& beam, const std::vector<size_t>& indices, const std::vector<Vec3>& coordinates,
        const std::vector<double>& depths, const std::vector<double>& distances, const GridInterpolant2D* kernels,
        const std::vector<double>& weights, const std::array<double, 3>& betas, double mu, double sad,
        double cutoff, bool divergent);
    std::vector<double> result();
private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace optirad