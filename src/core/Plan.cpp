#include "Plan.hpp"

namespace optirad {

void Plan::setName(const std::string& name) { m_name = name; }
const std::string& Plan::getName() const { return m_name; }

void Plan::addBeam(const Beam& beam) { m_beams.push_back(beam); }
const std::vector<Beam>& Plan::getBeams() const { return m_beams; }
size_t Plan::getNumBeams() const { return m_beams.size(); }

} // namespace optirad
