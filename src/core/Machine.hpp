#pragma once

#include <string>

namespace optirad {

class Machine {
public:
    void setName(const std::string& name);
    const std::string& getName() const;

    void setSAD(double sad);
    double getSAD() const;

private:
    std::string m_name;
    double m_sad = 1000.0; // Source-to-axis distance in mm
};

} // namespace optirad
