#include "Beam.hpp"

namespace optirad {

void Beam::setGantryAngle(double angle) { m_gantryAngle = angle; }
void Beam::setCouchAngle(double angle) { m_couchAngle = angle; }
void Beam::setIsocenter(double x, double y, double z) {
    m_isocenter[0] = x;
    m_isocenter[1] = y;
    m_isocenter[2] = z;
}

double Beam::getGantryAngle() const { return m_gantryAngle; }
double Beam::getCouchAngle() const { return m_couchAngle; }

} // namespace optirad
