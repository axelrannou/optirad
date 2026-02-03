#pragma once

#include "Structure.hpp"
#include <vector>
#include <optional>

namespace optirad {

class StructureSet {
public:
    void addStructure(const Structure& structure);
    size_t getNumStructures() const;

    const Structure& getStructure(size_t index) const;
    std::optional<const Structure*> findByName(const std::string& name) const;

private:
    std::vector<Structure> m_structures;
};

} // namespace optirad
