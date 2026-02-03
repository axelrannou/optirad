#include "DicomExporter.hpp"
#include "utils/Logger.hpp"

namespace optirad {

bool DicomExporter::exportRTPlan(const Plan& plan, const std::string& outputPath) {
    Logger::warn("DICOM RT Plan export not yet implemented");
    return false;
}

bool DicomExporter::exportRTDose(const DoseMatrix& dose, const std::string& outputPath) {
    Logger::warn("DICOM RT Dose export not yet implemented");
    return false;
}

} // namespace optirad
