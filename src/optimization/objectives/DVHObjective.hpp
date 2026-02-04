#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

/**
 * DVH-based objective: penalizes dose above/below threshold for a volume fraction
 * Example: D95 >= 95% means 95% of volume should receive at least prescription dose
 */
class DVHObjective : public ObjectiveFunction {
public:
    enum class Type {
        MIN_DOSE,  // At least volumeFraction should receive >= doseThreshold
        MAX_DOSE   // At most volumeFraction should receive >= doseThreshold
    };
    
    std::string getName() const override;
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
    
    void setType(Type type) { m_type = type; }
    void setDoseThreshold(double dose) { m_doseThreshold = dose; }
    void setVolumeFraction(double fraction) { m_volumeFraction = fraction; }
    
private:
    Type m_type = Type::MIN_DOSE;
    double m_doseThreshold = 0.0;
    double m_volumeFraction = 0.95;  // e.g., D95
};

} // namespace optirad
