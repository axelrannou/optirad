#pragma once

#include <string>

namespace optirad {

/// DICOM context captured during import.
/// Stores the UIDs needed to cross-reference the exported RT Plan with the
/// original study (CT, RT Struct) so that a compatible DICOM dataset is produced.
struct DicomContext {
    // Study-level UIDs (from CT slices)
    std::string studyInstanceUID;
    std::string frameOfReferenceUID;

    // RT Struct SOP Instance UID — to populate the Referenced Structure Set Sequence
    std::string rtStructSOPInstanceUID;

    // Study date/time (from CT, YYYYMMDD / HHMMSS)
    std::string studyDate;
    std::string studyTime;

    // Patient demographic tags (from CT)
    std::string patientBirthDate;
    std::string patientSex;

    /// Returns true when the minimum required UIDs are present.
    bool isValid() const { return !studyInstanceUID.empty(); }
};

} // namespace optirad
