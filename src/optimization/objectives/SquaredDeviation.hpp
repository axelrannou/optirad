#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

// (d - d_prescribed)^2 for target
class SquaredDeviation : public ObjectiveFunction {
public:
    std::string getName() const override;
    
    void setPrescribedDose(double dose);
    bool gpuDescription(DoseObjectiveSpec& output) const override {
        output = {0, m_prescribedDose, m_weight, 0.0, getActiveIndices()};
        return true;
    }

    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;

private:
    double m_prescribedDose = 0.0;
};

} // namespace optirad
