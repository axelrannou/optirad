#pragma once

#include "IDataExporter.hpp"
#include "dose/DoseMatrix.hpp"
#include <string>
#include <vector>

namespace optirad {

/// DICOM RT Plan exporter using DCMTK.
///
/// Writes a valid RT Plan DICOM file (SOP Class 1.2.840.10008.5.1.4.1.1.481.5)
/// from a completed leaf sequencing pipeline result. The exported file:
///   - Carries over the StudyInstanceUID and FrameOfReferenceUID from the original import,
///     so the RT Plan is part of the same study as the CT and RT Struct.
///   - Contains one Beam per entry in the Stf with full gantry/couch geometry.
///   - Encodes each aperture segment as a pair of control points (step-and-shoot).
///   - Stores MLC leaf positions in bankA/bankB order matching the DICOM MLCX convention.
///
/// Requires DCMTK at build time (OPTIRAD_HAS_DCMTK). Falls back to a no-op stub
/// that logs an error when DCMTK is not available.
class DicomExporter : public IDataExporter {
public:
    DicomExporter() = default;

    bool exportRTPlan(const Plan& plan,
                      const Stf& stf,
                      const std::vector<LeafSequenceResult>& sequences,
                      const DicomContext& context,
                      const std::string& outputDir,
                      std::string* outSOPInstanceUID = nullptr) override;

    /// Export a dose volume as a DICOM RT Dose file.
    /// @param dose       Dose matrix (values in Gy).
    /// @param context    DICOM context from the original import (for UIDs).
    /// @param outputDir  Directory where the RD.*.dcm file will be written.
    /// @param label      Optional comment string stored in DoseComment tag.
    /// @param referencedRTPlanSOPInstanceUID  Optional SOP Instance UID of an RT Plan
    ///                                        exported alongside this dose; when set,
    ///                                        a Referenced RT Plan Sequence is written
    ///                                        (required by the RT Dose IOD when
    ///                                        DoseSummationType is PLAN).
    bool exportRTDose(const DoseMatrix& dose,
                      const DicomContext& context,
                      const std::string& outputDir,
                      const std::string& label = "",
                      const std::string& referencedRTPlanSOPInstanceUID = "");

    /// Generate a new DICOM UID using DCMTK's UID generator.
    /// Returns an empty string when DCMTK is not available.
    static std::string generateUID();
};

} // namespace optirad
