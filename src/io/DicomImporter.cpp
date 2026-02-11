#include "DicomImporter.hpp"
#include "RTStructParser.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"
#include "utils/Logger.hpp"
#include <algorithm>

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
        
        // Rasterize contours to voxel indices
        if (patientData->getCTVolume()) {
            structures->rasterizeContours(patientData->getCTVolume()->getGrid());
        }
        
        patientData->setStructureSet(std::move(structures));
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
    grid.setDimensions(static_cast<size_t>(cols), static_cast<size_t>(rows), m_ctFiles.size());
    grid.setSpacing(spacingX, spacingY, spacingZ);
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
    
    RTStructParser parser;
    return parser.parse(m_rtStructFile.string());
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

} // namespace optirad
