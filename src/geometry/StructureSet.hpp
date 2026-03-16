#pragma once

#include "Structure.hpp"
#include <vector>
#include <memory>
#include <string>

namespace optirad {

class Grid; // Forward declaration

class StructureSet {
public:
    StructureSet() = default;
    ~StructureSet() = default;
    
    // Move-only (unique_ptr is not copyable)
    StructureSet(const StructureSet&) = delete;
    StructureSet& operator=(const StructureSet&) = delete;
    StructureSet(StructureSet&&) = default;
    StructureSet& operator=(StructureSet&&) = default;
    
    void addStructure(std::unique_ptr<Structure> structure);
    const Structure* getStructure(size_t index) const;
    Structure* getStructure(size_t index);
    const Structure* getStructureByName(const std::string& name) const;
    size_t getCount() const { return m_structures.size(); }
    
    /// Check if any structure has the given type (e.g., "EXTERNAL", "PTV")
    bool hasStructureOfType(const std::string& type) const {
        for (const auto& s : m_structures) {
            if (s && s->getType() == type) return true;
        }
        return false;
    }
    
    void clear() { m_structures.clear(); }

    /// Rasterize all structure contours to voxel indices
    void rasterizeContours(const Grid& ctGrid);

private:
    std::vector<std::unique_ptr<Structure>> m_structures;
};

} // namespace optirad
