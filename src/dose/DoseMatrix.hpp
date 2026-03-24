#pragma once

#include "Grid.hpp"
#include <vector>
#include <cstddef>

namespace optirad {

class DoseMatrix {
public:
    void setGrid(const Grid& grid);
    void allocate();

    double& at(size_t i, size_t j, size_t k);
    const double& at(size_t i, size_t j, size_t k) const;

    double getMax() const;
    double getMean() const;

    /// Trilinear interpolation at fractional voxel coordinates in the dose grid.
    /// Returns 0 if the point is outside the grid.
    double interpolateAt(double fi, double fj, double fk) const;

    const Grid& getGrid() const { return m_grid; }

    double* data();
    const double* data() const;
    size_t size() const;

private:
    Grid m_grid;
    std::vector<double> m_data;
};

} // namespace optirad
