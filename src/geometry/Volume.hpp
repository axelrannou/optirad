#pragma once

#include "Grid.hpp"
#include <vector>
#include <cstdint>

namespace optirad {

template<typename T>
class Volume {
public:
    void setGrid(const Grid& grid);
    const Grid& getGrid() const;

    void allocate();

    T& at(size_t i, size_t j, size_t k);
    const T& at(size_t i, size_t j, size_t k) const;

    T* data();
    const T* data() const;

    size_t size() const { return m_data.size(); }

private:
    Grid m_grid;
    std::vector<T> m_data;
};

// Typedefs for common volume types
using CTVolume = Volume<int16_t>;     // CT in HU
using DoseVolume = Volume<double>;    // Dose in Gy

} // namespace optirad
