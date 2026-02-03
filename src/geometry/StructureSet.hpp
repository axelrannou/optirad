#pragma once

#include "Structure.hpp"
#include <vector>
#include <memory>
#include <string>

namespace optirad {

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
    
    void clear() { m_structures.clear(); }

private:
    std::vector<std::unique_ptr<Structure>> m_structures;
};

} // namespace optirad
