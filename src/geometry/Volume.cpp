#include "Volume.hpp"
#include <stdexcept>

namespace optirad {

template<typename T>
void Volume<T>::setGrid(const Grid& grid) { m_grid = grid; }

template<typename T>
const Grid& Volume<T>::getGrid() const { return m_grid; }

template<typename T>
void Volume<T>::allocate() {
    m_data.resize(m_grid.getNumVoxels());
}

template<typename T>
T& Volume<T>::at(size_t i, size_t j, size_t k) {
    auto dims = m_grid.getDimensions();
    // Bounds checking
    if (i >= dims[0] || j >= dims[1] || k >= dims[2]) {
        throw std::out_of_range("Volume::at - index out of bounds");
    }
    return m_data[i + j * dims[0] + k * dims[0] * dims[1]];
}

template<typename T>
const T& Volume<T>::at(size_t i, size_t j, size_t k) const {
    auto dims = m_grid.getDimensions();
    // Bounds checking
    if (i >= dims[0] || j >= dims[1] || k >= dims[2]) {
        throw std::out_of_range("Volume::at - index out of bounds");
    }
    return m_data[i + j * dims[0] + k * dims[0] * dims[1]];
}

template<typename T>
T* Volume<T>::data() { return m_data.data(); }

template<typename T>
const T* Volume<T>::data() const { return m_data.data(); }

// Explicit instantiations
template class Volume<int16_t>;
template class Volume<double>;

} // namespace optirad
