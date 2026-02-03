#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

// max(0, d_min - d)^2 for target
class SquaredUnderdose : public ObjectiveFunction {
public:
    std::string getName() const override;
    
    void setMinDose(double dose);

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

private:
    double m_minDose = 0.0;
};

} // namespace optirad
