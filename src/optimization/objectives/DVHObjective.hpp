#pragma once

#include "../ObjectiveFunction.hpp"

namespace optirad {

/**
 * DVH-based objective: penalizes dose above/below threshold for a volume fraction.
 * 
 * MinDVH: At least volumeFraction of the structure should receive >= doseThreshold.
 *         Example: D95 >= 60 Gy → volumeFraction=0.95, doseThreshold=60
 * 
 * MaxDVH: At most volumeFraction of the structure should receive >= doseThreshold.
 *         Example: V20Gy <= 30% → volumeFraction=0.30, doseThreshold=20
 *         Example: Dmax <= 68 Gy → volumeFraction=0.02, doseThreshold=68
 */
class DVHObjective : public ObjectiveFunction {
public:
    enum class Type {
        MIN_DVH,  // At least volumeFraction should receive >= doseThreshold
        MAX_DVH   // At most volumeFraction should receive >= doseThreshold
    };
    
    std::string getName() const override;
    double compute(const std::vector<double>& dose) const override;
    std::vector<double> gradient(const std::vector<double>& dose) const override;
    
    void setType(Type type) { m_type = type; }
    Type getType() const { return m_type; }
    void setDoseThreshold(double dose) { m_doseThreshold = dose; }
    double getDoseThreshold() const { return m_doseThreshold; }
    void setVolumeFraction(double fraction) { m_volumeFraction = fraction; }
    double getVolumeFraction() const { return m_volumeFraction; }
    
private:
    Type m_type = Type::MIN_DVH;
    double m_doseThreshold = 0.0;
    double m_volumeFraction = 0.95;  // e.g., D95
};

} // namespace optirad
