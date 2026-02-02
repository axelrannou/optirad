#pragma once

#include <string>

namespace optirad {

enum class ConstraintType {
    MinDose,
    MaxDose,
    MeanDose
};

class Constraint {
public:
    void setType(ConstraintType type);
    void setValue(double value);
    void setStructureName(const std::string& name);

    ConstraintType getType() const;
    double getValue() const;
    const std::string& getStructureName() const;

    bool isSatisfied(double actualValue) const;

private:
    ConstraintType m_type = ConstraintType::MaxDose;
    double m_value = 0.0;
    std::string m_structureName;
};

} // namespace optirad
