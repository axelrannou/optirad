#include "Machine.hpp"

namespace optirad {

void Machine::setName(const std::string& name) { m_name = name; }
const std::string& Machine::getName() const { return m_name; }
void Machine::setSAD(double sad) { m_sad = sad; }
double Machine::getSAD() const { return m_sad; }

} // namespace optirad
