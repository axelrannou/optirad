#include "ObjectiveFunction.hpp"

namespace optirad {

void ObjectiveFunction::setWeight(double weight) { m_weight = weight; }
double ObjectiveFunction::getWeight() const { return m_weight; }
void ObjectiveFunction::setStructure(const Structure* structure) { m_structure = structure; }

} // namespace optirad
