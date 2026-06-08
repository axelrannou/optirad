#pragma once

#include "DicomContext.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/Aperture.hpp"
#include <string>
#include <vector>

namespace optirad {

/// Abstract interface for exporting treatment planning data to a DICOM file.
class IDataExporter {
public:
    virtual ~IDataExporter() = default;

    /// Export an RT Plan DICOM file to the given output directory.
    /// @param plan         The treatment plan (machine, fractions, beam geometry).
    /// @param stf          The steering file containing per-beam geometry.
    /// @param sequences    Leaf sequencing results (one per beam, ordered by beamIndex).
    /// @param context      DICOM UIDs captured during the original import.
    /// @param outputDir    Directory to write the RT Plan DICOM file into.
    /// @return             true on success, false on failure.
    virtual bool exportRTPlan(const Plan& plan,
                              const Stf& stf,
                              const std::vector<LeafSequenceResult>& sequences,
                              const DicomContext& context,
                              const std::string& outputDir) = 0;
};

} // namespace optirad
