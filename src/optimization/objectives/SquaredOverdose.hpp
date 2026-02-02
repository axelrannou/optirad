#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

// max(0, d - d_max)^2 for OAR
class SquaredOverdose : public ObjectiveFunction {
public:
    std::string getName() const override;
    
    void setMaxDose(double dose);

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

private:
    double m_maxDose = 0.0;
};

} // namespace optirad
