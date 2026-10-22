#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

// max(0, d_min - d)^2 for target
class SquaredUnderdose : public ObjectiveFunction {
public:
    std::string getName() const override;
    
    void setMinDose(double dose);
    bool gpuDescription(DoseObjectiveSpec& output) const override {
        output = {2, m_minDose, m_weight, 0.0, getActiveIndices()};
        return true;
    }

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

private:
    double m_minDose = 0.0;
};

} // namespace optirad
