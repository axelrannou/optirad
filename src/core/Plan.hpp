#pragma once

#include "Beam.hpp"
#include <string>
#include <vector>

namespace optirad {

class Plan {
public:
    void setName(const std::string& name);
    const std::string& getName() const;

    void addBeam(const Beam& beam);
    const std::vector<Beam>& getBeams() const;
    size_t getNumBeams() const;

private:
    std::string m_name;
    std::vector<Beam> m_beams;
};

} // namespace optirad
