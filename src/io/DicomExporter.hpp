#pragma once

#include "IDataExporter.hpp"

namespace optirad {

class DicomExporter : public IDataExporter {
public:
    bool exportDose(const DoseVolume& dose, const std::string& path) override;
    bool exportPlan(const Plan& plan, const std::string& path) override;
};

} // namespace optirad
