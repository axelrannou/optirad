#include "DicomImporter.hpp"
#include "Logger.hpp"

namespace optirad {

bool DicomImporter::load(const std::string& path) {
    Logger::info("Loading DICOM from: " + path);
    // TODO: Implement DICOM loading with DCMTK
    return true;
}

Patient DicomImporter::getPatient() const { return m_patient; }
CTVolume DicomImporter::getCT() const { return m_ct; }
StructureSet DicomImporter::getStructureSet() const { return m_structureSet; }

} // namespace optirad
