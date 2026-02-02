#pragma once

#include "Plan.hpp"
#include "Grid.hpp"
#include "DoseMatrix.hpp"
#include "DoseInfluenceMatrix.hpp"
#include <string>

namespace optirad {

class IDoseEngine {
public:
    virtual ~IDoseEngine() = default;

    virtual std::string getName() const = 0;

    // Calculate final dose distribution
    virtual DoseMatrix calculateDose(const Plan& plan, const Grid& grid) = 0;

    // Calculate dose influence matrix (Dij) for optimization
    virtual DoseInfluenceMatrix calculateDij(const Plan& plan, const Grid& grid) = 0;
};

} // namespace optirad
