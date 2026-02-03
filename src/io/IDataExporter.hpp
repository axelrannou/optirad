#pragma once

#include "Plan.hpp"
#include "Volume.hpp"
#include <string>

namespace optirad {

class IDataExporter {
public:
    virtual ~IDataExporter() = default;

    virtual bool exportDose(const DoseVolume& dose, const std::string& path) = 0;
    virtual bool exportPlan(const Plan& plan, const std::string& path) = 0;
};

} // namespace optirad
