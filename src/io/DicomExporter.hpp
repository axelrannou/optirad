#pragma once

#include <string>

namespace optirad {

class Plan;
class DoseMatrix;

class DicomExporter {
public:
    bool exportRTPlan(const Plan& plan, const std::string& outputPath);
    bool exportRTDose(const DoseMatrix& dose, const std::string& outputPath);
};

} // namespace optirad
