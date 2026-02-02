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

    double* data();
    const double* data() const;
    size_t size() const;

private:
    Grid m_grid;
    std::vector<double> m_data;
};

} // namespace optirad
