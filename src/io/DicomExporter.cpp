#include "DicomExporter.hpp"
#include "utils/Logger.hpp"

#include <filesystem>
#include <sstream>
#include <iomanip>
#include <ctime>
#include <algorithm>
#include <cstring>

#ifdef OPTIRAD_HAS_DCMTK
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmdata/dcuid.h>
#include <dcmtk/dcmdata/dcfilefo.h>
#include <dcmtk/dcmdata/dcdeftag.h>
#include <dcmtk/dcmdata/dcsequen.h>
#include <dcmtk/dcmdata/dcitem.h>
#include <dcmtk/dcmdata/dcvrus.h>
#endif

namespace optirad {

// ─────────────────────── helpers ───────────────────────

/// Format a double as a DICOM DS (Decimal String) value.
static std::string toDS(double v, int precision = 6) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(precision) << v;
    // Strip trailing zeros after decimal point
    std::string s = ss.str();
    auto dot = s.find('.');
    if (dot != std::string::npos) {
        size_t last = s.find_last_not_of('0');
        if (last != std::string::npos && last > dot)
            s = s.substr(0, last + 1);
        else if (last == dot)
            s = s.substr(0, dot);
    }
    return s;
}

/// Format a vector of doubles as a backslash-separated DICOM DS multi-value.
static std::string toDSMulti(const std::vector<double>& vals, int precision = 6) {
    std::string out;
    for (size_t i = 0; i < vals.size(); ++i) {
        if (i > 0) out += '\\';
        out += toDS(vals[i], precision);
    }
    return out;
}

/// YYYYMMDD from current local time.
static std::string currentDate() {
    std::time_t t = std::time(nullptr);
    std::tm* tm = std::localtime(&t);
    char buf[16];
    std::strftime(buf, sizeof(buf), "%Y%m%d", tm);
    return buf;
}

/// HHMMSS from current local time.
static std::string currentTime() {
    std::time_t t = std::time(nullptr);
    std::tm* tm = std::localtime(&t);
    char buf[16];
    std::strftime(buf, sizeof(buf), "%H%M%S", tm);
    return buf;
}

// ─────────────────────── public API ───────────────────────

std::string DicomExporter::generateUID() {
#ifdef OPTIRAD_HAS_DCMTK
    char uid[128];
    dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT);
    return std::string(uid);
#else
    return "";
#endif
}

bool DicomExporter::exportRTPlan(const Plan& plan,
                                  const Stf& stf,
                                  const std::vector<LeafSequenceResult>& sequences,
                                  const DicomContext& context,
                                  const std::string& outputDir)
{
#ifndef OPTIRAD_HAS_DCMTK
    Logger::error("DicomExporter: DCMTK not available — cannot export RT Plan");
    return false;
#else
    namespace fs = std::filesystem;

    if (sequences.empty()) {
        Logger::error("DicomExporter: no leaf sequencing results provided");
        return false;
    }
    if (stf.isEmpty()) {
        Logger::error("DicomExporter: Stf is empty");
        return false;
    }

    // Ensure output directory exists
    if (!fs::exists(outputDir)) {
        std::error_code ec;
        if (!fs::create_directories(outputDir, ec)) {
            Logger::error("DicomExporter: cannot create output directory: " + outputDir
                          + " (" + ec.message() + ")");
            return false;
        }
    }

    // ── Build a map from beamIndex → LeafSequenceResult ──
    std::vector<const LeafSequenceResult*> seqByBeam(stf.getCount(), nullptr);
    for (const auto& seq : sequences) {
        if (seq.beamIndex < stf.getCount())
            seqByBeam[seq.beamIndex] = &seq;
    }

    // ── SOP Instance UID (new for this file) ──
    const std::string sopInstanceUID = generateUID();
    const std::string seriesInstanceUID = generateUID();
    const std::string today = currentDate();
    const std::string now   = currentTime();

    // ── Patient info ──
    std::string patientName, patientID;
    if (plan.getPatientData() && plan.getPatientData()->getPatient()) {
        patientName = plan.getPatientData()->getPatient()->getName();
        patientID   = plan.getPatientData()->getPatient()->getID();
    }

    // ─────────────────────── DICOM dataset ───────────────────────
    DcmFileFormat fileFormat;
    DcmDataset* ds = fileFormat.getDataset();

    // ── Patient Module ──
    ds->putAndInsertString(DCM_PatientName,      patientName.c_str());
    ds->putAndInsertString(DCM_PatientID,        patientID.c_str());
    ds->putAndInsertString(DCM_PatientBirthDate, context.patientBirthDate.c_str());
    ds->putAndInsertString(DCM_PatientSex,       context.patientSex.c_str());

    // ── General Study Module ──
    std::string studyUID = context.studyInstanceUID.empty() ? generateUID() : context.studyInstanceUID;
    ds->putAndInsertString(DCM_StudyInstanceUID, studyUID.c_str());
    ds->putAndInsertString(DCM_StudyDate,        context.studyDate.empty() ? today.c_str() : context.studyDate.c_str());
    ds->putAndInsertString(DCM_StudyTime,        context.studyTime.empty() ? now.c_str()   : context.studyTime.c_str());
    ds->putAndInsertString(DCM_AccessionNumber,  "");
    ds->putAndInsertString(DCM_StudyID,          "1");

    // ── RT Series Module ──
    ds->putAndInsertString(DCM_SeriesInstanceUID, seriesInstanceUID.c_str());
    ds->putAndInsertString(DCM_SeriesNumber,      "1");
    ds->putAndInsertString(DCM_Modality,          "RTPLAN");

    // ── General Equipment Module ──
    ds->putAndInsertString(DCM_Manufacturer, "OptiRad");

    // ── Frame of Reference Module ──
    std::string forUID = context.frameOfReferenceUID.empty() ? generateUID() : context.frameOfReferenceUID;
    ds->putAndInsertString(DCM_FrameOfReferenceUID,       forUID.c_str());
    ds->putAndInsertString(DCM_PositionReferenceIndicator, "");

    // ── SOP Common Module ──
    ds->putAndInsertString(DCM_SOPClassUID,           UID_RTPlanStorage);
    ds->putAndInsertString(DCM_SOPInstanceUID,        sopInstanceUID.c_str());
    ds->putAndInsertString(DCM_InstanceCreationDate,  today.c_str());
    ds->putAndInsertString(DCM_InstanceCreationTime,  now.c_str());
    ds->putAndInsertString(DCM_SpecificCharacterSet,  "ISO_IR 6");

    // ── RT General Plan Module ──
    ds->putAndInsertString(DCM_RTPlanLabel,    plan.getName().substr(0, 16).c_str()); // max 16 chars
    ds->putAndInsertString(DCM_RTPlanName,     plan.getName().c_str());
    ds->putAndInsertString(DCM_RTPlanGeometry, "PATIENT");
    ds->putAndInsertString(DCM_RTPlanDate,     today.c_str());
    ds->putAndInsertString(DCM_RTPlanTime,     now.c_str());

    // ── Referenced Structure Set Sequence ──
    if (!context.rtStructSOPInstanceUID.empty()) {
        auto* refSSSeq = new DcmSequenceOfItems(DCM_ReferencedStructureSetSequence);
        auto* refSSItem = new DcmItem();
        refSSItem->putAndInsertString(DCM_ReferencedSOPClassUID,    UID_RTStructureSetStorage);
        refSSItem->putAndInsertString(DCM_ReferencedSOPInstanceUID, context.rtStructSOPInstanceUID.c_str());
        refSSSeq->insert(refSSItem);
        ds->insert(refSSSeq, OFTrue);
    }

    // ─────────────────────── Fraction Group Sequence ───────────────────────
    {
        auto* fgSeq = new DcmSequenceOfItems(DCM_FractionGroupSequence);
        auto* fgItem = new DcmItem();

        fgItem->putAndInsertString(DCM_FractionGroupNumber,        "1");
        fgItem->putAndInsertString(DCM_NumberOfFractionsPlanned,   std::to_string(plan.getNumOfFractions()).c_str());
        fgItem->putAndInsertString(DCM_NumberOfBeams,              std::to_string(stf.getCount()).c_str());
        fgItem->putAndInsertString(DCM_NumberOfBrachyApplicationSetups, "0");

        // Referenced Beam Sequence (one item per beam)
        auto* refBeamSeq = new DcmSequenceOfItems(DCM_ReferencedBeamSequence);
        for (size_t i = 0; i < stf.getCount(); ++i) {
            auto* rbItem = new DcmItem();
            rbItem->putAndInsertString(DCM_ReferencedBeamNumber,
                                       std::to_string(static_cast<int>(i) + 1).c_str());

            double meterset = 0.0;
            if (seqByBeam[i]) meterset = seqByBeam[i]->totalMU;
            rbItem->putAndInsertString(DCM_BeamMeterset, toDS(meterset).c_str());

            refBeamSeq->insert(rbItem);
        }
        fgItem->insert(refBeamSeq, OFTrue);

        fgSeq->insert(fgItem);
        ds->insert(fgSeq, OFTrue);
    }

    // ─────────────────────── Beam Sequence ───────────────────────
    {
        auto* beamSeq = new DcmSequenceOfItems(DCM_BeamSequence);

        for (size_t beamIdx = 0; beamIdx < stf.getCount(); ++beamIdx) {
            const Beam* beam = stf.getBeam(beamIdx);
            if (!beam) continue;

            const LeafSequenceResult* seq = seqByBeam[beamIdx];
            int numSegments = (seq && !seq->segments.empty())
                              ? static_cast<int>(seq->segments.size()) : 0;
            int numCP       = 2 * numSegments; // step-and-shoot: 2 CPs per segment

            auto* beamItem = new DcmItem();

            // ── Beam identification ──
            beamItem->putAndInsertString(DCM_BeamNumber,
                std::to_string(static_cast<int>(beamIdx) + 1).c_str());
            beamItem->putAndInsertString(DCM_BeamName,
                ("Beam_" + std::to_string(static_cast<int>(beamIdx) + 1)).c_str());
            std::string beamDesc = "Gantry " + toDS(beam->getGantryAngle(), 1) + " Couch " + toDS(beam->getCouchAngle(), 1);
            beamItem->putAndInsertString(DCM_BeamDescription, beamDesc.c_str());
            beamItem->putAndInsertString(DCM_BeamType,              "STATIC");
            beamItem->putAndInsertString(DCM_RadiationType,         "PHOTON");
            beamItem->putAndInsertString(DCM_TreatmentDeliveryType, "TREATMENT");
            beamItem->putAndInsertString(DCM_TreatmentMachineName,
                plan.getMachine().getName().substr(0, 16).c_str());
            beamItem->putAndInsertString(DCM_PrimaryDosimeterUnit,  "MU");
            beamItem->putAndInsertString(DCM_SourceAxisDistance,
                toDS(beam->getSAD()).c_str());

            // ── Counts (mandatory attributes) ──
            beamItem->putAndInsertString(DCM_NumberOfWedges,       "0");
            beamItem->putAndInsertString(DCM_NumberOfCompensators, "0");
            beamItem->putAndInsertString(DCM_NumberOfBoli,         "0");
            beamItem->putAndInsertString(DCM_NumberOfBlocks,       "0");
            beamItem->putAndInsertString(DCM_FinalCumulativeMetersetWeight,
                toDS(seq ? seq->totalMU : 0.0).c_str());
            beamItem->putAndInsertString(DCM_NumberOfControlPoints,
                std::to_string(numCP).c_str());

            // ── Beam Limiting Device Sequence (static geometry) ──
            {
                auto* bldSeq = new DcmSequenceOfItems(DCM_BeamLimitingDeviceSequence);

                // ASYMX jaw
                auto* asymx = new DcmItem();
                asymx->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMX");
                asymx->putAndInsertString(DCM_NumberOfLeafJawPairs,     "1");
                bldSeq->insert(asymx);

                // ASYMY jaw
                auto* asymy = new DcmItem();
                asymy->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMY");
                asymy->putAndInsertString(DCM_NumberOfLeafJawPairs,     "1");
                bldSeq->insert(asymy);

                // MLCX — include Leaf Position Boundaries if available
                if (seq && !seq->leafPairBoundariesZ.empty()) {
                    auto* mlcx = new DcmItem();
                    mlcx->putAndInsertString(DCM_RTBeamLimitingDeviceType, "MLCX");
                    int numLeafPairs = static_cast<int>(seq->leafPairBoundariesZ.size()) - 1;
                    mlcx->putAndInsertString(DCM_NumberOfLeafJawPairs,
                        std::to_string(numLeafPairs).c_str());
                    mlcx->putAndInsertString(DCM_LeafPositionBoundaries,
                        toDSMulti(seq->leafPairBoundariesZ).c_str());
                    bldSeq->insert(mlcx);
                } else if (seq && !seq->segments.empty() && !seq->segments[0].bankA.empty()) {
                    // Fallback: infer boundaries from bankA size and bixel width
                    int numLeafPairs = static_cast<int>(seq->segments[0].bankA.size());
                    auto* mlcx = new DcmItem();
                    mlcx->putAndInsertString(DCM_RTBeamLimitingDeviceType, "MLCX");
                    mlcx->putAndInsertString(DCM_NumberOfLeafJawPairs,
                        std::to_string(numLeafPairs).c_str());
                    bldSeq->insert(mlcx);
                }

                beamItem->insert(bldSeq, OFTrue);
            }

            // ── Control Point Sequence (2 CPs per aperture segment) ──
            if (numCP > 0 && seq) {
                auto* cpSeq = new DcmSequenceOfItems(DCM_ControlPointSequence);

                // Beam isocenter (from beam geometry)
                const Vec3& iso = beam->getIsocenter();
                std::string isoStr = toDS(iso[0]) + "\\" + toDS(iso[1]) + "\\" + toDS(iso[2]);

                // Nominal beam energy (from first ray)
                double energy = 6.0;
                if (!beam->getRays().empty())
                    energy = beam->getRays()[0].getEnergy();

                std::string gantryAngleStr = toDS(beam->getGantryAngle(), 4);
                std::string couchAngleStr  = toDS(beam->getCouchAngle(), 4);

                double cumulativeMU = 0.0;

                for (int s = 0; s < numSegments; ++s) {
                    const Aperture& seg = seq->segments[static_cast<size_t>(s)];
                    double segMU = seg.weight;

                    // Compute jaw extent from current segment's leaf positions
                    double jawXmin = 0.0, jawXmax = 0.0;
                    double jawYmin = 0.0, jawYmax = 0.0;

                    if (!seg.bankA.empty() && !seg.bankB.empty()) {
                        jawXmin = *std::min_element(seg.bankA.begin(), seg.bankA.end());
                        jawXmax = *std::max_element(seg.bankB.begin(), seg.bankB.end());
                    }
                    if (!seq->leafPairBoundariesZ.empty()) {
                        jawYmin = seq->leafPairBoundariesZ.front();
                        jawYmax = seq->leafPairBoundariesZ.back();
                    }

                    // Build MLCX position string: [bankA[0..N-1], bankB[0..N-1]]
                    std::string mlcPositions;
                    {
                        std::vector<double> allPos;
                        allPos.insert(allPos.end(), seg.bankA.begin(), seg.bankA.end());
                        allPos.insert(allPos.end(), seg.bankB.begin(), seg.bankB.end());
                        mlcPositions = toDSMulti(allPos, 2);
                    }

                    // ── CP start (index 2*s) ──
                    {
                        auto* cp = new DcmItem();
                        cp->putAndInsertString(DCM_ControlPointIndex,
                            std::to_string(2 * s).c_str());
                        // Full beam geometry is only required in the first CP,
                        // but we write it in every start CP for robustness.
                        cp->putAndInsertString(DCM_NominalBeamEnergy,     toDS(energy, 2).c_str());
                        cp->putAndInsertString(DCM_GantryAngle,           gantryAngleStr.c_str());
                        cp->putAndInsertString(DCM_GantryRotationDirection,"NONE");
                        cp->putAndInsertString(DCM_BeamLimitingDeviceAngle,"0");
                        cp->putAndInsertString(DCM_BeamLimitingDeviceRotationDirection, "NONE");
                        cp->putAndInsertString(DCM_PatientSupportAngle,   couchAngleStr.c_str());
                        cp->putAndInsertString(DCM_PatientSupportRotationDirection, "NONE");
                        cp->putAndInsertString(DCM_TableTopEccentricAngle,"0");
                        cp->putAndInsertString(DCM_TableTopEccentricRotationDirection, "NONE");
                        cp->putAndInsertString(DCM_TableTopVerticalPosition,    "");
                        cp->putAndInsertString(DCM_TableTopLongitudinalPosition, "");
                        cp->putAndInsertString(DCM_TableTopLateralPosition,     "");
                        cp->putAndInsertString(DCM_IsocenterPosition,     isoStr.c_str());
                        cp->putAndInsertString(DCM_CumulativeMetersetWeight,
                            toDS(cumulativeMU).c_str());

                        // Beam Limiting Device Position Sequence
                        auto* bldpSeq = new DcmSequenceOfItems(DCM_BeamLimitingDevicePositionSequence);

                        auto* jawX = new DcmItem();
                        jawX->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMX");
                        jawX->putAndInsertString(DCM_LeafJawPositions,
                            (toDS(jawXmin, 2) + "\\" + toDS(jawXmax, 2)).c_str());
                        bldpSeq->insert(jawX);

                        auto* jawY = new DcmItem();
                        jawY->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMY");
                        jawY->putAndInsertString(DCM_LeafJawPositions,
                            (toDS(jawYmin, 2) + "\\" + toDS(jawYmax, 2)).c_str());
                        bldpSeq->insert(jawY);

                        if (!seg.bankA.empty()) {
                            auto* mlcx = new DcmItem();
                            mlcx->putAndInsertString(DCM_RTBeamLimitingDeviceType, "MLCX");
                            mlcx->putAndInsertString(DCM_LeafJawPositions, mlcPositions.c_str());
                            bldpSeq->insert(mlcx);
                        }

                        cp->insert(bldpSeq, OFTrue);
                        cpSeq->insert(cp);
                    }

                    cumulativeMU += segMU;

                    // ── CP end (index 2*s + 1) ──
                    {
                        auto* cp = new DcmItem();
                        cp->putAndInsertString(DCM_ControlPointIndex,
                            std::to_string(2 * s + 1).c_str());
                        cp->putAndInsertString(DCM_CumulativeMetersetWeight,
                            toDS(cumulativeMU).c_str());

                        // Repeat device positions at end of segment (required by some viewers)
                        auto* bldpSeq = new DcmSequenceOfItems(DCM_BeamLimitingDevicePositionSequence);

                        auto* jawX = new DcmItem();
                        jawX->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMX");
                        jawX->putAndInsertString(DCM_LeafJawPositions,
                            (toDS(jawXmin, 2) + "\\" + toDS(jawXmax, 2)).c_str());
                        bldpSeq->insert(jawX);

                        auto* jawY = new DcmItem();
                        jawY->putAndInsertString(DCM_RTBeamLimitingDeviceType, "ASYMY");
                        jawY->putAndInsertString(DCM_LeafJawPositions,
                            (toDS(jawYmin, 2) + "\\" + toDS(jawYmax, 2)).c_str());
                        bldpSeq->insert(jawY);

                        if (!seg.bankA.empty()) {
                            auto* mlcx = new DcmItem();
                            mlcx->putAndInsertString(DCM_RTBeamLimitingDeviceType, "MLCX");
                            mlcx->putAndInsertString(DCM_LeafJawPositions, mlcPositions.c_str());
                            bldpSeq->insert(mlcx);
                        }

                        cp->insert(bldpSeq, OFTrue);
                        cpSeq->insert(cp);
                    }
                }

                beamItem->insert(cpSeq, OFTrue);
            }

            beamSeq->insert(beamItem);
        }

        ds->insert(beamSeq, OFTrue);
    }

    // ── Write the file ──
    std::string filename = (fs::path(outputDir) / ("RP." + sopInstanceUID + ".dcm")).string();

    OFCondition status = fileFormat.saveFile(filename.c_str(), EXS_LittleEndianExplicit,
                                             EET_ExplicitLength, EGL_recalcGL,
                                             EPD_withoutPadding,
                                             OFstatic_cast(Uint32, 0),
                                             OFstatic_cast(Uint32, 0),
                                             EWM_fileformat);
    if (status.bad()) {
        Logger::error("DicomExporter: failed to save RT Plan: " + std::string(status.text()));
        return false;
    }

    Logger::info("DicomExporter: RT Plan written to " + filename);
    return true;
#endif
}

bool DicomExporter::exportRTDose(const DoseMatrix& dose,
                                  const DicomContext& context,
                                  const std::string& outputDir,
                                  const std::string& label)
{
#ifndef OPTIRAD_HAS_DCMTK
    Logger::error("DicomExporter: DCMTK not available — cannot export RT Dose");
    return false;
#else
    namespace fs = std::filesystem;

    const Grid& grid = dose.getGrid();
    auto dims = grid.getDimensions();
    const size_t nx = dims[0], ny = dims[1], nz = dims[2];
    if (nx == 0 || ny == 0 || nz == 0) {
        Logger::error("DicomExporter: empty dose grid");
        return false;
    }

    // Ensure output directory exists
    if (!fs::exists(outputDir)) {
        std::error_code ec;
        if (!fs::create_directories(outputDir, ec)) {
            Logger::error("DicomExporter: cannot create output dir: " + outputDir
                          + " (" + ec.message() + ")");
            return false;
        }
    }

    const std::string sopInstanceUID  = generateUID();
    const std::string seriesInstanceUID = generateUID();
    const std::string today = currentDate();
    const std::string now   = currentTime();

    // ── Dose grid scaling: scale max dose to near UINT32_MAX ──
    double maxDose = dose.getMax();
    if (maxDose <= 0.0) maxDose = 1.0;
    const double doseGridScaling = maxDose / static_cast<double>(0xFFFFFFFEu);

    // ── Build uint32 pixel array (frame-major, row-major, column-minor) ──
    std::vector<Uint32> pixelData(nx * ny * nz);
    for (size_t k = 0; k < nz; ++k) {
        for (size_t j = 0; j < ny; ++j) {
            for (size_t i = 0; i < nx; ++i) {
                double d = dose.at(i, j, k);
                pixelData[k * ny * nx + j * nx + i] =
                    static_cast<Uint32>(d / doseGridScaling + 0.5);
            }
        }
    }

    // ── GridFrameOffsetVector (relative z-offset of each frame from first) ──
    auto zCoords = grid.getZCoordinates();
    const double z0 = zCoords.empty() ? 0.0 : zCoords[0];
    const double dz = grid.getSpacing()[2];
    std::vector<double> frameOffsets(nz);
    for (size_t k = 0; k < nz; ++k)
        frameOffsets[k] = (k < zCoords.size() ? zCoords[k] : z0 + k * dz) - z0;

    // ────────────────────────── Dataset ──────────────────────────
    DcmFileFormat fileFormat;
    DcmDataset* ds = fileFormat.getDataset();

    // Patient Module
    ds->putAndInsertString(DCM_PatientName,      "");
    ds->putAndInsertString(DCM_PatientID,        "");
    ds->putAndInsertString(DCM_PatientBirthDate, context.patientBirthDate.c_str());
    ds->putAndInsertString(DCM_PatientSex,       context.patientSex.c_str());

    // General Study Module
    std::string studyUID = context.studyInstanceUID.empty() ? generateUID() : context.studyInstanceUID;
    ds->putAndInsertString(DCM_StudyInstanceUID, studyUID.c_str());
    ds->putAndInsertString(DCM_StudyDate,
        context.studyDate.empty() ? today.c_str() : context.studyDate.c_str());
    ds->putAndInsertString(DCM_StudyTime,
        context.studyTime.empty() ? now.c_str()   : context.studyTime.c_str());
    ds->putAndInsertString(DCM_AccessionNumber,  "");
    ds->putAndInsertString(DCM_StudyID,          "1");

    // RT Series
    ds->putAndInsertString(DCM_SeriesInstanceUID, seriesInstanceUID.c_str());
    ds->putAndInsertString(DCM_SeriesNumber,      "2");
    ds->putAndInsertString(DCM_Modality,          "RTDOSE");

    // General Equipment
    ds->putAndInsertString(DCM_Manufacturer, "OptiRad");

    // Frame of Reference
    std::string forUID = context.frameOfReferenceUID.empty() ? generateUID() : context.frameOfReferenceUID;
    ds->putAndInsertString(DCM_FrameOfReferenceUID,        forUID.c_str());
    ds->putAndInsertString(DCM_PositionReferenceIndicator, "");

    // SOP Common
    ds->putAndInsertString(DCM_SOPClassUID,           UID_RTDoseStorage);
    ds->putAndInsertString(DCM_SOPInstanceUID,        sopInstanceUID.c_str());
    ds->putAndInsertString(DCM_InstanceCreationDate,  today.c_str());
    ds->putAndInsertString(DCM_InstanceCreationTime,  now.c_str());
    ds->putAndInsertString(DCM_SpecificCharacterSet,  "ISO_IR 6");

    // Image Pixel Module
    ds->putAndInsertString(DCM_SamplesPerPixel,          "1");
    ds->putAndInsertString(DCM_PhotometricInterpretation, "MONOCHROME2");
    ds->putAndInsertString(DCM_NumberOfFrames,           std::to_string(nz).c_str());
    ds->putAndInsertUint16(DCM_Rows,                     static_cast<Uint16>(ny));
    ds->putAndInsertUint16(DCM_Columns,                  static_cast<Uint16>(nx));
    ds->putAndInsertUint16(DCM_BitsAllocated,            32);
    ds->putAndInsertUint16(DCM_BitsStored,               32);
    ds->putAndInsertUint16(DCM_HighBit,                  31);
    ds->putAndInsertUint16(DCM_PixelRepresentation,      0);

    // Image geometry
    auto spacing = grid.getSpacing();
    // PixelSpacing = row spacing (y-direction), col spacing (x-direction)
    ds->putAndInsertString(DCM_PixelSpacing,
        (toDS(spacing[1], 6) + "\\" + toDS(spacing[0], 6)).c_str());

    auto origin = grid.getOrigin();
    ds->putAndInsertString(DCM_ImagePositionPatient,
        (toDS(origin[0]) + "\\" + toDS(origin[1]) + "\\" + toDS(origin[2])).c_str());

    const auto& orient = grid.getImageOrientation();
    ds->putAndInsertString(DCM_ImageOrientationPatient,
        (toDS(orient[0]) + "\\" + toDS(orient[1]) + "\\" + toDS(orient[2]) + "\\" +
         toDS(orient[3]) + "\\" + toDS(orient[4]) + "\\" + toDS(orient[5])).c_str());

    ds->putAndInsertString(DCM_SliceThickness, toDS(spacing[2]).c_str());

    // Frame Increment Pointer → GridFrameOffsetVector (tag 3004,000C)
    {
        auto* fip = new DcmAttributeTag(DCM_FrameIncrementPointer);
        fip->putTagVal(DCM_GridFrameOffsetVector, 0);
        ds->insert(fip, OFTrue);
    }
    ds->putAndInsertString(DCM_GridFrameOffsetVector, toDSMulti(frameOffsets, 4).c_str());

    // RT Dose Module
    ds->putAndInsertString(DCM_DoseUnits,         "GY");
    ds->putAndInsertString(DCM_DoseType,          "PHYSICAL");
    ds->putAndInsertString(DCM_DoseSummationType, "PLAN");
    ds->putAndInsertString(DCM_DoseGridScaling,   toDS(doseGridScaling, 10).c_str());
    if (!label.empty())
        ds->putAndInsertString(DCM_DoseComment, label.c_str());

    // Pixel data (uint32 values)
    ds->putAndInsertUint32Array(DCM_PixelData, pixelData.data(),
                                static_cast<unsigned long>(pixelData.size()));

    // ── Write the file ──
    std::string filename = (fs::path(outputDir) / ("RD." + sopInstanceUID + ".dcm")).string();
    OFCondition status = fileFormat.saveFile(filename.c_str(), EXS_LittleEndianExplicit,
                                             EET_ExplicitLength, EGL_recalcGL,
                                             EPD_withoutPadding,
                                             OFstatic_cast(Uint32, 0),
                                             OFstatic_cast(Uint32, 0),
                                             EWM_fileformat);
    if (status.bad()) {
        Logger::error("DicomExporter: failed to save RT Dose: " + std::string(status.text()));
        return false;
    }

    Logger::info("DicomExporter: RT Dose written to " + filename);
    return true;
#endif
}

} // namespace optirad
