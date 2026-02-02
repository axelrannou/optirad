#pragma once

#include "Structure.hpp"
#include <vector>
#include <string>

namespace optirad {

class ObjectiveFunction {
public:
    virtual ~ObjectiveFunction() = default;

    virtual std::string getName() const = 0;

    // Compute objective value given dose distribution
    virtual double compute(const std::vector<double>& dose) const = 0;

    // Compute gradient with respect to dose
    virtual std::vector<double> gradient(const std::vector<double>& dose) const = 0;

    void setWeight(double weight);
    double getWeight() const;

    void setStructure(const Structure* structure);

protected:
    double m_weight = 1.0;
    const Structure* m_structure = nullptr;
};

} // namespace optirad
