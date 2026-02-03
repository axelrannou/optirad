#include "DicomExporter.hpp"
#include "Logger.hpp"

namespace optirad {

bool DicomExporter::exportDose(const DoseVolume& dose, const std::string& path) {
    Logger::info("Exporting DICOM RT Dose to: " + path);
    // TODO: Implement DICOM RT Dose export
    return true;
}

bool DicomExporter::exportPlan(const Plan& plan, const std::string& path) {
    Logger::info("Exporting DICOM RT Plan to: " + path);
    // TODO: Implement DICOM RT Plan export
    return true;
}

} // namespace optirad
