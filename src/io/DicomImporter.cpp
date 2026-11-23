#include "DicomImporter.hpp"
#include "DicomContext.hpp"
#include "RTStructParser.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"
#include "dose/DoseMatrix.hpp"
#include "segmentation/BodyContourGenerator.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <map>
#include <optional>
#include <sstream>

#ifdef OPTIRAD_HAS_DCMTK
#include <dcmtk/dcmdata/dctk.h>
#include <dcmtk/dcmimgle/dcmimage.h>
#endif

namespace optirad {

// DICOM SOP Class UIDs
namespace SopClass {
    constexpr const char* CTImageStorage = "1.2.840.10008.5.1.4.1.1.2";
    constexpr const char* RTStructureSet = "1.2.840.10008.5.1.4.1.1.481.3";
    constexpr const char* RTPlan = "1.2.840.10008.5.1.4.1.1.481.5";
    constexpr const char* RTDose = "1.2.840.10008.5.1.4.1.1.481.2";
}

bool DicomImporter::canImport(const std::string& path) const {
    namespace fs = std::filesystem;
    if (!fs::exists(path)) {
        Logger::warn("Path does not exist: " + path);
        return false;
    }
    
    if (fs::is_directory(path)) {
        int fileCount = 0;
        for (const auto& entry : fs::recursive_directory_iterator(path)) {
            if (entry.is_regular_file()) {
                auto filename = entry.path().filename().string();
                // Check length before substr to prevent buffer underflow
                if (filename.length() >= 2 && filename.substr(0, 2) == "._") {
                    continue;
                }
                fileCount++;
            }
        }
        Logger::info("Found " + std::to_string(fileCount) + " files in directory");
        return fileCount > 0;
    }
    
    return fs::is_regular_file(path);
}

std::unique_ptr<PatientData> DicomImporter::importAll(const std::string& dirPath) {
    if (!loadDirectory(dirPath)) {
        Logger::error("Failed to load DICOM directory: " + dirPath);
        return nullptr;
    }
    
    auto patientData = std::make_unique<PatientData>();
    
    // Import patient info
    patientData->setPatient(importPatient(dirPath));
    
    // Import CT volume
    auto ctVolume = importCTVolume();
    if (ctVolume) {
        Logger::info("CT volume loaded: " + std::to_string(ctVolume->size()) + " voxels");
        patientData->setCTVolume(std::move(ctVolume));
        
        // Convert to electron density
        patientData->convertHUtoED();
        Logger::info("Converted HU to electron density");
    }
    
    // Import structures with contours
    auto structures = importStructuresWithContours();
    if (structures) {
        Logger::info("Loaded " + std::to_string(structures->getCount()) + " structures");
        
        // Auto-generate BODY contour if no EXTERNAL structure exists
        if (!structures->hasStructureOfType("EXTERNAL") && patientData->getCTVolume()) {
            Logger::info("No BODY/EXTERNAL structure found, auto-generating from CT threshold");
            auto bodyStructure = BodyContourGenerator::generate(*patientData->getCTVolume());
            if (bodyStructure) {
                structures->addStructure(std::move(bodyStructure));
            }
        }
        
        // Rasterize contours to voxel indices
        if (patientData->getCTVolume()) {
            structures->rasterizeContours(patientData->getCTVolume()->getGrid());
        }
        
        patientData->setStructureSet(std::move(structures));
    }

    // Capture RT Struct SOPInstanceUID for the DICOM context (used during export)
#ifdef OPTIRAD_HAS_DCMTK
    if (!m_rtStructFile.empty()) {
        DcmFileFormat rtStructFF;
        if (rtStructFF.loadFile(m_rtStructFile.c_str()).good()) {
            OFString rtStructSopUID;
            if (rtStructFF.getDataset()->findAndGetOFString(DCM_SOPInstanceUID, rtStructSopUID).good()) {
                m_context.rtStructSOPInstanceUID = rtStructSopUID.c_str();
            }
        }
    }
#endif

    // Import RT Dose if available — cache it for the caller to take via takeImportedDose()
    auto [doseMatrix, doseGrid] = importRTDose();
    if (doseMatrix) {
        m_importedDoseMatrix = doseMatrix;
        m_importedDoseGrid   = doseGrid;
        Logger::info("RT Dose imported: max=" + std::to_string(doseMatrix->getMax()) + " Gy");
        // Grid spacing is stored as (row, column, slice) = (y, x, z)
        const auto sp = doseGrid->getSpacing();
        patientData->setImportedDoseGridResolution({sp[1], sp[0], sp[2]});
    }

    // Import RT Plan if available. A plan that cannot be read never blocks the rest of the import.
    if (!m_rtPlanFile.empty()) {
        if (auto rtPlan = importRTPlan()) {
            for (const auto& w : rtPlan->warnings) Logger::warn("RTPLAN: " + w);
            Logger::info("RT Plan imported: '" + rtPlan->label + "', " + std::to_string(rtPlan->beams.size()) +
                         " treatment beam(s), " + std::to_string(rtPlan->numFractions) + " fraction(s)");
            patientData->setImportedPlan(std::move(rtPlan));
        }
    }
    
    return patientData;
}

std::unique_ptr<Volume<int16_t>> DicomImporter::importCTVolume() {
#ifdef OPTIRAD_HAS_DCMTK
    if (m_ctFiles.empty()) {
        Logger::warn("No CT files to import");
        return nullptr;
    }
    
    sortCTFilesByPosition();
    
    // Read first slice to get dimensions and geometric info
    DcmFileFormat firstFile;
    if (!firstFile.loadFile(m_ctFiles[0].c_str()).good()) {
        Logger::error("Failed to load first CT slice");
        return nullptr;
    }
    
    DcmDataset* ds = firstFile.getDataset();
    if (!ds) {
        Logger::error("Failed to get dataset from first CT slice");
        return nullptr;
    }
    
    Uint16 rows = 0, cols = 0;
    ds->findAndGetUint16(DCM_Rows, rows);
    ds->findAndGetUint16(DCM_Columns, cols);
    
    // Get pixel spacing
    OFString pixelSpacing;
    double spacingX = 1.0, spacingY = 1.0;
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 0).good()) {
        try {
            spacingY = std::stod(pixelSpacing.c_str());
        } catch (const std::exception& e) {
            Logger::warn("Invalid pixel spacing value (Y): " + std::string(pixelSpacing.c_str()));
            spacingY = 1.0;
        }
    }
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 1).good()) {
        try {
            spacingX = std::stod(pixelSpacing.c_str());
        } catch (const std::exception& e) {
            Logger::warn("Invalid pixel spacing value (X): " + std::string(pixelSpacing.c_str()));
            spacingX = 1.0;
        }
    }
    
    // Get image position (origin)
    double originX = 0, originY = 0, originZ = 0;
    OFString pos;
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 0).good()) {
        try {
            originX = std::stod(pos.c_str());
        } catch (const std::exception& e) {
            Logger::warn("Invalid image position (X): " + std::string(pos.c_str()));
            originX = 0.0;
        }
    }
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 1).good()) {
        try {
            originY = std::stod(pos.c_str());
        } catch (const std::exception& e) {
            Logger::warn("Invalid image position (Y): " + std::string(pos.c_str()));
            originY = 0.0;
        }
    }
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 2).good()) {
        try {
            originZ = std::stod(pos.c_str());
        } catch (const std::exception& e) {
            Logger::warn("Invalid image position (Z): " + std::string(pos.c_str()));
            originZ = 0.0;
        }
    }
    
    // Get patient position (HFS, FFS, HFP, etc.)
    OFString patientPosition;
    ds->findAndGetOFString(DCM_PatientPosition, patientPosition);
    std::string patientPosStr = patientPosition.c_str();
    if (patientPosStr.empty()) {
        patientPosStr = "HFS";  // Default: Head First Supine
    }
    
    // Get image orientation (patient) - 6 values: row direction + column direction
    std::array<double, 6> imageOrientation = {1, 0, 0, 0, 1, 0};  // Default: standard axial
    OFString orientStr;
    for (int i = 0; i < 6; ++i) {
        if (ds->findAndGetOFString(DCM_ImageOrientationPatient, orientStr, i).good()) {
            try {
                imageOrientation[i] = std::stod(orientStr.c_str());
            } catch (const std::exception& e) {
                Logger::warn("Invalid image orientation value at index " + std::to_string(i));
                // Keep default value
            }
        }
    }
    
    // Get slice thickness (may differ from calculated spacing)
    double sliceThickness = 0.0;
    Float64 thickness = 0.0;
    if (ds->findAndGetFloat64(DCM_SliceThickness, thickness).good()) {
        sliceThickness = thickness;
    }
    
    // Calculate actual slice spacing from first two slices
    double spacingZ = 1.0;
    if (m_ctFiles.size() > 1) {
        DcmFileFormat secondFile;
        if (secondFile.loadFile(m_ctFiles[1].c_str()).good()) {
            double z2 = 0;
            if (secondFile.getDataset()->findAndGetOFString(DCM_ImagePositionPatient, pos, 2).good()) {
                z2 = std::stod(pos.c_str());
            }
            spacingZ = std::abs(z2 - originZ);
        }
    }
    
    // If slice thickness not found, use calculated spacing
    if (sliceThickness == 0.0) {
        sliceThickness = spacingZ;
    }
    
    // Create grid with all geometric information
    Grid grid;
    grid.setDimensions(static_cast<size_t>(rows), static_cast<size_t>(cols), m_ctFiles.size());
    grid.setSpacing(spacingY, spacingX, spacingZ);
    grid.setOrigin({originX, originY, originZ});
    grid.setPatientPosition(patientPosStr);
    grid.setImageOrientation(imageOrientation);
    grid.setSliceThickness(sliceThickness);
    
    Logger::info("CT Grid: " + std::to_string(cols) + "x" + std::to_string(rows) + "x" + std::to_string(m_ctFiles.size()));
    Logger::info("Spacing: " + std::to_string(spacingX) + "x" + std::to_string(spacingY) + "x" + std::to_string(spacingZ) + " mm");
    Logger::info("Patient Position: " + patientPosStr);
    Logger::info("Image Orientation: [" + 
                std::to_string(imageOrientation[0]) + ", " + std::to_string(imageOrientation[1]) + ", " + std::to_string(imageOrientation[2]) + ", " +
                std::to_string(imageOrientation[3]) + ", " + std::to_string(imageOrientation[4]) + ", " + std::to_string(imageOrientation[5]) + "]");
    Logger::info("Slice Thickness: " + std::to_string(sliceThickness) + " mm");
    
    // Create volume
    auto volume = std::make_unique<Volume<int16_t>>();
    volume->setGrid(grid);
    volume->allocate();
    
    // Load all slices
    size_t sliceSize = static_cast<size_t>(rows) * cols;
    for (size_t z = 0; z < m_ctFiles.size(); ++z) {
        DcmFileFormat dcmFile;
        if (!dcmFile.loadFile(m_ctFiles[z].c_str()).good()) {
            Logger::warn("Failed to load slice " + std::to_string(z));
            continue;
        }
        
        DcmDataset* sliceDs = dcmFile.getDataset();
        if (!sliceDs) {
            Logger::warn("Failed to get dataset for slice " + std::to_string(z));
            continue;
        }
        
        // Get rescale slope and intercept
        double rescaleSlope = 1.0, rescaleIntercept = 0.0;
        if (!sliceDs->findAndGetFloat64(DCM_RescaleSlope, rescaleSlope).good()) {
            Logger::warn("Missing RescaleSlope for slice " + std::to_string(z) + ", using default 1.0");
        }
        if (!sliceDs->findAndGetFloat64(DCM_RescaleIntercept, rescaleIntercept).good()) {
            Logger::warn("Missing RescaleIntercept for slice " + std::to_string(z) + ", using default 0.0");
        }
        
        // Get pixel data
        const Uint16* pixelData = nullptr;
        unsigned long count = 0;
        if (sliceDs->findAndGetUint16Array(DCM_PixelData, pixelData, &count).good() && pixelData) {
            int16_t* destPtr = volume->data() + z * sliceSize;
            for (size_t i = 0; i < sliceSize && i < count; ++i) {
                // Apply rescale to get HU
                double hu = static_cast<double>(pixelData[i]) * rescaleSlope + rescaleIntercept;
                destPtr[i] = static_cast<int16_t>(std::clamp(hu, -32768.0, 32767.0));
            }
        }
    }
    
    return volume;
#else
    Logger::warn("DCMTK not available - cannot import CT volume");
    return nullptr;
#endif
}

std::unique_ptr<StructureSet> DicomImporter::importStructuresWithContours() {
    if (m_rtStructFile.empty()) {
        Logger::warn("No RT Structure Set file found");
        return nullptr;
    }
    
    // RTStructParser now uses TBB internally for parallel processing
    RTStructParser parser;
    auto structureSet = parser.parse(m_rtStructFile.string());
    
    if (structureSet) {
        Logger::info("Successfully parsed structure set with " + 
                    std::to_string(structureSet->getCount()) + " structures");
    }
    
    return structureSet;
}

void DicomImporter::sortCTFilesByPosition() {
#ifdef OPTIRAD_HAS_DCMTK
    std::vector<std::pair<double, std::filesystem::path>> filesWithZ;
    
    for (const auto& file : m_ctFiles) {
        DcmFileFormat dcmFile;
        if (dcmFile.loadFile(file.c_str()).good()) {
            OFString pos;
            double z = 0;
            if (dcmFile.getDataset()->findAndGetOFString(DCM_ImagePositionPatient, pos, 2).good()) {
                z = std::stod(pos.c_str());
            }
            filesWithZ.emplace_back(z, file);
        }
    }
    
    std::sort(filesWithZ.begin(), filesWithZ.end());
    
    m_ctFiles.clear();
    for (const auto& [z, path] : filesWithZ) {
        m_ctFiles.push_back(path);
    }
#endif
}

std::unique_ptr<Patient> DicomImporter::importPatient(const std::string& path) {
    if (m_ctFiles.empty() && !loadDirectory(path)) {
        Logger::error("Failed to load DICOM directory: " + path);
        return nullptr;
    }
    
    auto patient = std::make_unique<Patient>();
    
#ifdef OPTIRAD_HAS_DCMTK
    if (!m_ctFiles.empty()) {
        DcmFileFormat fileFormat;
        if (fileFormat.loadFile(m_ctFiles[0].c_str()).good()) {
            DcmDataset* dataset = fileFormat.getDataset();

            OFString patientName, patientID;
            dataset->findAndGetOFString(DCM_PatientName, patientName);
            dataset->findAndGetOFString(DCM_PatientID, patientID);

            patient->setName(patientName.c_str());
            patient->setID(patientID.c_str());
            Logger::info("Loaded patient: " + std::string(patientName.c_str()));

            // Capture DICOM context UIDs for later export
            OFString studyUID, forUID, studyDate, studyTime, birthDate, sex;
            dataset->findAndGetOFString(DCM_StudyInstanceUID,    studyUID);
            dataset->findAndGetOFString(DCM_FrameOfReferenceUID, forUID);
            dataset->findAndGetOFString(DCM_StudyDate,           studyDate);
            dataset->findAndGetOFString(DCM_StudyTime,           studyTime);
            dataset->findAndGetOFString(DCM_PatientBirthDate,    birthDate);
            dataset->findAndGetOFString(DCM_PatientSex,          sex);

            m_context.studyInstanceUID  = studyUID.c_str();
            m_context.frameOfReferenceUID = forUID.c_str();
            m_context.studyDate         = studyDate.c_str();
            m_context.studyTime         = studyTime.c_str();
            m_context.patientBirthDate  = birthDate.c_str();
            m_context.patientSex        = sex.c_str();
        }
    }
#else
    Logger::warn("DCMTK not available - using stub patient data");
    patient->setName("StubPatient");
    patient->setID("STUB001");
#endif
    
    return patient;
}

std::unique_ptr<StructureSet> DicomImporter::importStructures(const std::string& path) {
    if (m_rtStructFile.empty()) {
        loadDirectory(path);
    }
    return importStructuresWithContours();
}

bool DicomImporter::loadDirectory(const std::string& dirPath) {
    namespace fs = std::filesystem;
    
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
        Logger::error("Invalid DICOM directory: " + dirPath);
        return false;
    }
    
    Logger::info("Scanning DICOM directory: " + dirPath);
    
    m_ctFiles.clear();
    m_rtStructFile.clear();
    m_rtPlanFile.clear();
    m_rtDoseFile.clear();
    
    scanDirectory(fs::path(dirPath));
    sortCTFilesByPosition();
    
    Logger::info("Scan complete:");
    Logger::info("  CT slices: " + std::to_string(m_ctFiles.size()));
    if (!m_rtStructFile.empty()) Logger::info("  RT Structure Set: " + m_rtStructFile.filename().string());
    if (!m_rtPlanFile.empty()) Logger::info("  RT Plan: " + m_rtPlanFile.filename().string());
    if (!m_rtDoseFile.empty()) Logger::info("  RT Dose: " + m_rtDoseFile.filename().string());
    
    return !m_ctFiles.empty();
}

void DicomImporter::scanDirectory(const std::filesystem::path& dir) {
    namespace fs = std::filesystem;
    
    for (const auto& entry : fs::recursive_directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        
        auto filename = entry.path().filename().string();
        if (filename[0] == '.' || filename.substr(0, 2) == "._" || filename == "DICOMDIR") continue;
        
        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        if (ext == ".txt" || ext == ".log" || ext == ".xml" || ext == ".json") continue;
        
        std::string sopClass = getSOPClassUID(entry.path());
        
        if (!sopClass.empty()) {
            if (sopClass == SopClass::CTImageStorage) {
                m_ctFiles.push_back(entry.path());
            } else if (sopClass == SopClass::RTStructureSet) {
                m_rtStructFile = entry.path();
            } else if (sopClass == SopClass::RTPlan) {
                m_rtPlanFile = entry.path();
            } else if (sopClass == SopClass::RTDose) {
                m_rtDoseFile = entry.path();
            }
        }
    }
}

std::string DicomImporter::getSOPClassUID(const std::filesystem::path& file) const {
#ifdef OPTIRAD_HAS_DCMTK
    DcmFileFormat fileFormat;
    if (fileFormat.loadFile(file.c_str()).good()) {
        OFString sopClass;
        if (fileFormat.getDataset()->findAndGetOFString(DCM_SOPClassUID, sopClass).good()) {
            return std::string(sopClass.c_str());
        }
    }
#endif
    return "";
}

bool DicomImporter::loadCTSeries(const std::string& dirPath) {
    return loadDirectory(dirPath) && !m_ctFiles.empty();
}

bool DicomImporter::loadRTStruct(const std::string& filePath) {
    m_rtStructFile = filePath;
    return std::filesystem::exists(filePath);
}

bool DicomImporter::loadRTPlan(const std::string& filePath) {
    m_rtPlanFile = filePath;
    return std::filesystem::exists(filePath);
}

bool DicomImporter::loadRTDose(const std::string& filePath) {
    m_rtDoseFile = filePath;
    return std::filesystem::exists(filePath);
}

namespace {

double wrap360(double a) {
    a = std::fmod(a, 360.0);
    return a < 0.0 ? a + 360.0 : a;
}

#ifdef OPTIRAD_HAS_DCMTK
std::vector<double> readDoubles(DcmItem* item, const DcmTagKey& tag) {
    std::vector<double> out;
    OFString str;
    if (!item || item->findAndGetOFStringArray(tag, str).bad()) return out;
    std::stringstream ss(str.c_str());
    std::string tok;
    while (std::getline(ss, tok, '\\')) {
        try { out.push_back(std::stod(tok)); } catch (...) {}
    }
    return out;
}

std::optional<double> readDouble(DcmItem* item, const DcmTagKey& tag) {
    const auto v = readDoubles(item, tag);
    if (v.empty()) return std::nullopt;
    return v.front();
}

std::string readString(DcmItem* item, const DcmTagKey& tag) {
    OFString s;
    if (!item || item->findAndGetOFString(tag, s).bad()) return {};
    std::string out = s.c_str();
    const auto first = out.find_first_not_of(' ');
    const auto last = out.find_last_not_of(' ');
    return first == std::string::npos ? std::string() : out.substr(first, last - first + 1);
}

std::vector<DcmItem*> sequenceItems(DcmItem* parent, const DcmTagKey& tag) {
    std::vector<DcmItem*> out;
    DcmSequenceOfItems* seq = nullptr;
    if (parent && parent->findAndGetSequence(tag, seq).good() && seq) {
        for (unsigned long i = 0; i < seq->card(); ++i) out.push_back(seq->getItem(i));
    }
    return out;
}
#endif

} // namespace

std::shared_ptr<ImportedPlan> DicomImporter::importRTPlan(const std::string& filePath) {
#ifdef OPTIRAD_HAS_DCMTK
    const std::filesystem::path path = filePath.empty() ? m_rtPlanFile : std::filesystem::path(filePath);
    if (path.empty()) return nullptr;

    DcmFileFormat fileFormat;
    if (!fileFormat.loadFile(path.c_str()).good()) {
        Logger::error("Failed to load RT Plan file: " + path.string());
        return nullptr;
    }
    DcmDataset* ds = fileFormat.getDataset();
    if (!ds) return nullptr;

    auto plan = std::make_shared<ImportedPlan>();
    plan->label = readString(ds, DCM_RTPlanLabel);
    plan->name = readString(ds, DCM_RTPlanName);
    plan->sopInstanceUID = readString(ds, DCM_SOPInstanceUID);

    const auto setups = sequenceItems(ds, DCM_PatientSetupSequence);
    if (!setups.empty()) {
        const std::string pos = readString(setups.front(), DCM_PatientPosition);
        if (!pos.empty()) plan->patientPosition = pos;
    }

    // Fraction group: number of fractions and the MU per fraction of each beam
    std::map<int, double> metersetByBeam;
    const auto fractionGroups = sequenceItems(ds, DCM_FractionGroupSequence);
    if (fractionGroups.size() > 1)
        plan->warnings.push_back("plan has " + std::to_string(fractionGroups.size()) +
                                 " fraction groups; only the first one is used");
    if (!fractionGroups.empty()) {
        plan->numFractions = static_cast<int>(readDouble(fractionGroups.front(), DCM_NumberOfFractionsPlanned).value_or(0.0));
        for (DcmItem* ref : sequenceItems(fractionGroups.front(), DCM_ReferencedBeamSequence)) {
            const auto num = readDouble(ref, DCM_ReferencedBeamNumber);
            const auto mu = readDouble(ref, DCM_BeamMeterset);
            if (num && mu) metersetByBeam[static_cast<int>(*num)] = *mu;
        }
    }

    // Prescription: first TARGET dose reference with a prescription dose, else any dose reference with one
    for (DcmItem* ref : sequenceItems(ds, DCM_DoseReferenceSequence)) {
        const auto dose = readDouble(ref, DCM_TargetPrescriptionDose);
        if (!dose || *dose <= 0.0) continue;
        const bool target = readString(ref, DCM_DoseReferenceType) == "TARGET";
        if (plan->prescribedDose <= 0.0 || target) plan->prescribedDose = *dose;
        if (target) break;
    }

    int skipped = 0;
    for (DcmItem* beamItem : sequenceItems(ds, DCM_BeamSequence)) {
        ImportedBeam beam;
        beam.number = static_cast<int>(readDouble(beamItem, DCM_BeamNumber).value_or(0.0));
        beam.name = readString(beamItem, DCM_BeamName);
        const std::string deliveryType = readString(beamItem, DCM_TreatmentDeliveryType);
        if (!deliveryType.empty() && deliveryType != "TREATMENT") { ++skipped; continue; }

        const auto meterset = metersetByBeam.find(beam.number);
        if (meterset == metersetByBeam.end() || meterset->second <= 0.0) {
            plan->warnings.push_back("beam '" + beam.name + "' has no beam meterset and is skipped");
            continue;
        }
        beam.metersetMU = meterset->second;
        beam.machineName = readString(beamItem, DCM_TreatmentMachineName);
        beam.radiationType = readString(beamItem, DCM_RadiationType);
        beam.sad = readDouble(beamItem, DCM_SourceAxisDistance).value_or(1000.0);

        for (DcmItem* dev : sequenceItems(beamItem, DCM_BeamLimitingDeviceSequence)) {
            const std::string type = readString(dev, DCM_RTBeamLimitingDeviceType);
            if (type == "MLCX") beam.leafBoundaries = readDoubles(dev, DCM_LeafPositionBoundaries);
            else if (type == "MLCY")
                plan->warnings.push_back("beam '" + beam.name + "' uses an MLCY collimator, which is not supported");
        }

        const auto fluence = sequenceItems(beamItem, DCM_PrimaryFluenceModeSequence);
        if (!fluence.empty() && readString(fluence.front(), DCM_FluenceMode) == "NON_STANDARD") {
            beam.fff = true;
            const std::string id = readString(fluence.front(), DCM_FluenceModeID);
            plan->warnings.push_back("beam '" + beam.name + "' uses a non-standard fluence mode" +
                                     (id.empty() ? std::string() : " (" + id + ")"));
        }

        // Control points: absent attributes keep the value of the previous control point
        ImportedControlPoint state;
        std::string direction;
        double collimator = 0.0, couch = 0.0;
        std::optional<double> initialCollimator, initialCouch;
        bool dynamicCollimatorOrCouch = false;
        const auto cpItems = sequenceItems(beamItem, DCM_ControlPointSequence);
        for (size_t k = 0; k < cpItems.size(); ++k) {
            DcmItem* cp = cpItems[k];
            if (auto g = readDouble(cp, DCM_GantryAngle)) state.gantryAngle = wrap360(*g);
            if (auto w = readDouble(cp, DCM_CumulativeMetersetWeight)) state.cumulativeWeight = *w;
            if (auto c = readDouble(cp, DCM_BeamLimitingDeviceAngle)) {
                if (k > 0 && std::abs(*c - collimator) > 1e-6) dynamicCollimatorOrCouch = true;
                if (!initialCollimator) initialCollimator = *c;
                collimator = *c;
            }
            if (auto c = readDouble(cp, DCM_PatientSupportAngle)) {
                if (k > 0 && std::abs(*c - couch) > 1e-6) dynamicCollimatorOrCouch = true;
                if (!initialCouch) initialCouch = *c;
                couch = *c;
            }
            if (k == 0) {
                direction = readString(cp, DCM_GantryRotationDirection);
                beam.nominalEnergy = readDouble(cp, DCM_NominalBeamEnergy).value_or(0.0);
                beam.doseRateMUPerMin = readDouble(cp, DCM_DoseRateSet).value_or(0.0);
                const auto iso = readDoubles(cp, DCM_IsocenterPosition);
                if (iso.size() == 3) beam.isocenter = {iso[0], iso[1], iso[2]};
                else plan->warnings.push_back("beam '" + beam.name + "' has no isocenter position");
            }
            for (DcmItem* dev : sequenceItems(cp, DCM_BeamLimitingDevicePositionSequence)) {
                const std::string type = readString(dev, DCM_RTBeamLimitingDeviceType);
                const auto pos = readDoubles(dev, DCM_LeafJawPositions);
                if ((type == "X" || type == "ASYMX") && pos.size() == 2) state.jawX = {pos[0], pos[1]};
                else if ((type == "Y" || type == "ASYMY") && pos.size() == 2) state.jawY = {pos[0], pos[1]};
                else if (type == "MLCX" && pos.size() >= 2 && pos.size() % 2 == 0) {
                    const size_t n = pos.size() / 2;
                    state.bankA.assign(pos.begin(), pos.begin() + static_cast<long>(n));
                    state.bankB.assign(pos.begin() + static_cast<long>(n), pos.end());
                }
            }
            beam.controlPoints.push_back(state);
        }
        beam.collimatorAngle = initialCollimator.value_or(collimator);
        couch = wrap360(initialCouch.value_or(couch));
        beam.couchAngle = couch > 180.0 ? couch - 360.0 : couch;
        if (dynamicCollimatorOrCouch)
            plan->warnings.push_back("beam '" + beam.name + "' rotates the collimator or couch during delivery; "
                                     "only the first control point's angles are used");

        if (beam.controlPoints.empty()) {
            plan->warnings.push_back("beam '" + beam.name + "' has no control points and is skipped");
            continue;
        }
        const double finalWeight = readDouble(beamItem, DCM_FinalCumulativeMetersetWeight)
                                       .value_or(beam.controlPoints.back().cumulativeWeight);
        if (finalWeight <= 0.0) {
            plan->warnings.push_back("beam '" + beam.name + "' has a zero final meterset weight and is skipped");
            continue;
        }
        for (auto& cp : beam.controlPoints) cp.cumulativeWeight /= finalWeight;

        beam.isArc = (direction == "CW" || direction == "CC") && beam.controlPoints.size() >= 2;
        beam.clockwise = (direction != "CC");
        if (!beam.controlPoints.front().bankA.empty() && beam.leafBoundaries.size() != beam.controlPoints.front().bankA.size() + 1)
            plan->warnings.push_back("beam '" + beam.name + "' has leaf positions that do not match its leaf boundaries");

        plan->beams.push_back(std::move(beam));
    }

    if (skipped > 0)
        Logger::info("RTPLAN: skipped " + std::to_string(skipped) + " setup/non-treatment beam(s)");
    if (plan->beams.empty()) {
        Logger::warn("RT Plan has no treatment beams: " + path.string());
        return nullptr;
    }
    return plan;
#else
    (void)filePath;
    Logger::warn("DCMTK not available: RT Plan import disabled");
    return nullptr;
#endif
}

std::pair<std::shared_ptr<DoseMatrix>, std::shared_ptr<Grid>> DicomImporter::importRTDose() {
#ifdef OPTIRAD_HAS_DCMTK
    if (m_rtDoseFile.empty()) {
        return {nullptr, nullptr};
    }

    DcmFileFormat fileFormat;
    if (!fileFormat.loadFile(m_rtDoseFile.c_str()).good()) {
        Logger::error("Failed to load RT Dose file: " + m_rtDoseFile.string());
        return {nullptr, nullptr};
    }

    DcmDataset* ds = fileFormat.getDataset();
    if (!ds) {
        Logger::error("RT Dose file has no dataset");
        return {nullptr, nullptr};
    }

    // Get dimensions
    Uint16 rows = 0, cols = 0;
    ds->findAndGetUint16(DCM_Rows, rows);
    ds->findAndGetUint16(DCM_Columns, cols);

    // Number of frames (slices)
    OFString numFramesStr;
    long numFrames = 1;
    if (ds->findAndGetOFString(DCM_NumberOfFrames, numFramesStr).good()) {
        numFrames = std::atol(numFramesStr.c_str());
    }
    if (numFrames < 1) numFrames = 1;

    Logger::info("RT Dose dimensions: " + std::to_string(cols) + "x" +
                std::to_string(rows) + "x" + std::to_string(numFrames));

    // Pixel spacing
    double spacingX = 1.0, spacingY = 1.0;
    OFString pixelSpacing;
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 0).good()) {
        try { spacingY = std::stod(pixelSpacing.c_str()); } catch (...) {}
    }
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 1).good()) {
        try { spacingX = std::stod(pixelSpacing.c_str()); } catch (...) {}
    }

    // Image position (origin)
    double originX = 0, originY = 0, originZ = 0;
    OFString pos;
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 0).good()) {
        try { originX = std::stod(pos.c_str()); } catch (...) {}
    }
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 1).good()) {
        try { originY = std::stod(pos.c_str()); } catch (...) {}
    }
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 2).good()) {
        try { originZ = std::stod(pos.c_str()); } catch (...) {}
    }

    // Grid frame offset vector → compute slice spacing
    // Read offset[0] and offset[1] to get the inter-frame distance.
    // Offsets are relative to ImagePositionPatient of the first frame.
    double spacingZ = 1.0;
    bool framesDescending = false;
    if (numFrames > 1) {
        double off0 = 0.0, off1 = 0.0;
        OFString offsetStr;
        if (ds->findAndGetOFString(DCM_GridFrameOffsetVector, offsetStr, 0).good()) {
            try { off0 = std::stod(offsetStr.c_str()); } catch (...) {}
        }
        if (ds->findAndGetOFString(DCM_GridFrameOffsetVector, offsetStr, 1).good()) {
            try { off1 = std::stod(offsetStr.c_str()); } catch (...) {}
        }
        double step = off1 - off0;
        if (std::abs(step) > 0.001) {
            framesDescending = (step < 0);
            spacingZ = std::abs(step);
        }
    }

    // Dose grid scaling factor
    Float64 doseGridScaling = 1.0;
    ds->findAndGetFloat64(DCM_DoseGridScaling, doseGridScaling);

    // BitsAllocated determines pixel word size (typically 32 for RT Dose)
    Uint16 bitsAllocated = 32;
    ds->findAndGetUint16(DCM_BitsAllocated, bitsAllocated);

    Logger::info("RT Dose: DoseGridScaling=" + std::to_string(doseGridScaling) +
                ", BitsAllocated=" + std::to_string(bitsAllocated));

    // Image orientation (defaults to standard axial if absent)
    std::array<double, 6> imageOrientation = {1, 0, 0, 0, 1, 0};
    OFString orientStr;
    for (int i = 0; i < 6; ++i) {
        if (ds->findAndGetOFString(DCM_ImageOrientationPatient, orientStr, i).good()) {
            try { imageOrientation[i] = std::stod(orientStr.c_str()); } catch (...) {}
        }
    }

    // Create grid — if frames are descending, adjust origin to the last frame
    // so the grid always goes in the positive slice direction.
    auto grid = std::make_shared<Grid>();
    grid->setDimensions(static_cast<size_t>(rows), static_cast<size_t>(cols),
                        static_cast<size_t>(numFrames));
    grid->setSpacing(spacingY, spacingX, spacingZ);

    if (framesDescending && numFrames > 1) {
        // Origin becomes position of last DICOM frame (which is lowest Z)
        double adjustedOriginZ = originZ + (numFrames - 1) * (-spacingZ); // offset from first frame
        grid->setOrigin({originX, originY, adjustedOriginZ});
    } else {
        grid->setOrigin({originX, originY, originZ});
    }
    grid->setImageOrientation(imageOrientation);

    Logger::info("RT Dose grid: origin=(" + std::to_string(grid->getOrigin()[0]) + ", " +
                std::to_string(grid->getOrigin()[1]) + ", " + std::to_string(grid->getOrigin()[2]) +
                "), spacing=(" + std::to_string(spacingX) + ", " + std::to_string(spacingY) +
                ", " + std::to_string(spacingZ) + ")" +
                (framesDescending ? " [frames reversed]" : ""));

    size_t totalVoxels = static_cast<size_t>(rows) * cols * numFrames;
    size_t sliceSize = static_cast<size_t>(rows) * cols;

    auto doseMatrix = std::make_shared<DoseMatrix>();
    doseMatrix->setGrid(*grid);
    doseMatrix->allocate();

    // Read pixel data via raw byte buffer to avoid DCMTK VR mismatch
    // (RT Dose pixel data is OW but contains 32-bit unsigned integers)
    DcmElement* pixelElement = nullptr;
    OFCondition findPxStatus = ds->findAndGetElement(DCM_PixelData, pixelElement);

    if (!findPxStatus.good() || !pixelElement) {
        Logger::error("RT Dose: pixel data element not found");
        return {nullptr, nullptr};
    }

    // DICOM stores pixels row-major: data[col + row*cols + frame*rows*cols]
    // DoseMatrix::at(i,j,k) uses column-major: data[i + j*rows + k*rows*cols]
    // We must transpose within each frame so at(row, col, frame) returns the
    // correct DICOM pixel at (row, col, frame).
    // Also handle descending frame order by reversing frame index.
    auto transposeAndStore = [&](auto* srcData) {
        for (long frame = 0; frame < numFrames; ++frame) {
            long dstFrame = framesDescending ? (numFrames - 1 - frame) : frame;
            for (size_t row = 0; row < rows; ++row) {
                for (size_t col = 0; col < cols; ++col) {
                    size_t dicomFlat = col + row * cols + frame * sliceSize;
                    // at(row, col, dstFrame) index:
                    size_t atFlat = row + col * rows + dstFrame * sliceSize;
                    doseMatrix->data()[atFlat] = static_cast<double>(srcData[dicomFlat]) * doseGridScaling;
                }
            }
        }
    };

    if (bitsAllocated == 32) {
        // Get raw byte buffer and reinterpret as 32-bit unsigned
        Uint8* rawBuffer = nullptr;
        pixelElement->getUint8Array(rawBuffer);
        if (rawBuffer) {
            const Uint32* data32 = reinterpret_cast<const Uint32*>(rawBuffer);
            transposeAndStore(data32);
        } else {
            Logger::error("RT Dose: failed to read 32-bit pixel data buffer");
            return {nullptr, nullptr};
        }
    } else {
        // 16-bit pixel data
        Uint16* data16 = nullptr;
        pixelElement->getUint16Array(data16);
        if (data16) {
            transposeAndStore(data16);
        } else {
            Logger::error("RT Dose: failed to read 16-bit pixel data");
            return {nullptr, nullptr};
        }
    }

    Logger::info("RT Dose loaded: " + std::to_string(totalVoxels) + " voxels, max=" +
                std::to_string(doseMatrix->getMax()) + " Gy");

    return {doseMatrix, grid};
#else
    Logger::warn("DCMTK not available - cannot import RT Dose");
    return {nullptr, nullptr};
#endif
}

} // namespace optirad
