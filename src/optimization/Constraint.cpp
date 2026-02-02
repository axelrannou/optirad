#include "Constraint.hpp"

namespace optirad {

void Constraint::setType(ConstraintType type) { m_type = type; }
void Constraint::setValue(double value) { m_value = value; }
void Constraint::setStructureName(const std::string& name) { m_structureName = name; }

ConstraintType Constraint::getType() const { return m_type; }
double Constraint::getValue() const { return m_value; }
const std::string& Constraint::getStructureName() const { return m_structureName; }

bool Constraint::isSatisfied(double actualValue) const {
    switch (m_type) {
        case ConstraintType::MinDose: return actualValue >= m_value;
        case ConstraintType::MaxDose: return actualValue <= m_value;
        case ConstraintType::MeanDose: return actualValue <= m_value;
        default: return false;
    }
}

} // namespace optirad
