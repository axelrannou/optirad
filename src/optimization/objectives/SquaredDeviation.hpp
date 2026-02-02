#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

// (d - d_prescribed)^2 for target
class SquaredDeviation : public ObjectiveFunction {
public:
    std::string getName() const override;
    
    void setPrescribedDose(double dose);

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

private:
    double m_prescribedDose = 0.0;
};

} // namespace optirad
