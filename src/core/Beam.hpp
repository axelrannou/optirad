#pragma once

#include <string>

namespace optirad {

class Beam {
public:
    void setGantryAngle(double angle);
    void setCouchAngle(double angle);
    void setIsocenter(double x, double y, double z);

    double getGantryAngle() const;
    double getCouchAngle() const;

private:
    double m_gantryAngle = 0.0;
    double m_couchAngle = 0.0;
    double m_isocenter[3] = {0.0, 0.0, 0.0};
};

} // namespace optirad
