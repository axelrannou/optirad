#include "Patient.hpp"

namespace optirad {

void Patient::setId(const std::string& id) { m_id = id; }
void Patient::setName(const std::string& name) { m_name = name; }
const std::string& Patient::getId() const { return m_id; }
const std::string& Patient::getName() const { return m_name; }

} // namespace optirad
