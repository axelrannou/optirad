#pragma once

#include "geometry/Grid.hpp"
#include <memory>
#include <vector>

namespace optirad {

class CudaRayTracer {
public:
    CudaRayTracer(const Grid& grid, const double* density);
    ~CudaRayTracer();
    std::vector<double> depths(const Vec3& source, const std::vector<Vec3>& targets);
private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace optirad