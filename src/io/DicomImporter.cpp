#include "DicomImporter.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"
#include "geometry/Grid.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <map>

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
                if (filename.substr(0, 2) == "._") continue;
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
    
    // Read first slice to get dimensions
    DcmFileFormat firstFile;
    if (!firstFile.loadFile(m_ctFiles[0].c_str()).good()) {
        Logger::error("Failed to load first CT slice");
        return nullptr;
    }
    
    DcmDataset* ds = firstFile.getDataset();
    
    Uint16 rows = 0, cols = 0;
    ds->findAndGetUint16(DCM_Rows, rows);
    ds->findAndGetUint16(DCM_Columns, cols);
    
    // Get pixel spacing
    OFString pixelSpacing;
    double spacingX = 1.0, spacingY = 1.0;
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 0).good()) {
        spacingY = std::stod(pixelSpacing.c_str());
    }
    if (ds->findAndGetOFString(DCM_PixelSpacing, pixelSpacing, 1).good()) {
        spacingX = std::stod(pixelSpacing.c_str());
    }
    
    // Get image position (origin)
    double originX = 0, originY = 0, originZ = 0;
    OFString pos;
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 0).good()) originX = std::stod(pos.c_str());
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 1).good()) originY = std::stod(pos.c_str());
    if (ds->findAndGetOFString(DCM_ImagePositionPatient, pos, 2).good()) originZ = std::stod(pos.c_str());
    
    // Calculate slice thickness from first two slices
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
    
    // Create grid - use 3 separate arguments as per Grid API
    Grid grid;
    grid.setDimensions(static_cast<size_t>(cols), static_cast<size_t>(rows), m_ctFiles.size());
    grid.setSpacing(spacingX, spacingY, spacingZ);
    grid.setOrigin({originX, originY, originZ});
    
    Logger::info("CT Grid: " + std::to_string(cols) + "x" + std::to_string(rows) + "x" + std::to_string(m_ctFiles.size()));
    Logger::info("Spacing: " + std::to_string(spacingX) + "x" + std::to_string(spacingY) + "x" + std::to_string(spacingZ) + " mm");
    
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
        
        // Get rescale slope and intercept
        double rescaleSlope = 1.0, rescaleIntercept = 0.0;
        sliceDs->findAndGetFloat64(DCM_RescaleSlope, rescaleSlope);
        sliceDs->findAndGetFloat64(DCM_RescaleIntercept, rescaleIntercept);
        
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
#ifdef OPTIRAD_HAS_DCMTK
    if (m_rtStructFile.empty()) {
        Logger::warn("No RT Structure Set file found");
        return nullptr;
    }
    
    auto structures = std::make_unique<StructureSet>();
    
    DcmFileFormat fileFormat;
    if (!fileFormat.loadFile(m_rtStructFile.c_str()).good()) {
        Logger::error("Failed to load RT Structure Set");
        return nullptr;
    }
    
    DcmDataset* dataset = fileFormat.getDataset();
    
    // Build ROI number to name mapping
    std::map<int, std::string> roiNames;
    DcmSequenceOfItems* roiSequence = nullptr;
    if (dataset->findAndGetSequence(DCM_StructureSetROISequence, roiSequence).good()) {
        for (unsigned long i = 0; i < roiSequence->card(); ++i) {
            DcmItem* item = roiSequence->getItem(i);
            if (!item) continue;
            
            OFString roiName;
            Sint32 roiNumber = 0;
            item->findAndGetOFString(DCM_ROIName, roiName);
            item->findAndGetSint32(DCM_ROINumber, roiNumber);
            roiNames[roiNumber] = roiName.c_str();
        }
    }
    
    // Read ROI contour data
    DcmSequenceOfItems* contourSequence = nullptr;
    if (dataset->findAndGetSequence(DCM_ROIContourSequence, contourSequence).good()) {
        for (unsigned long i = 0; i < contourSequence->card(); ++i) {
            DcmItem* roiItem = contourSequence->getItem(i);
            if (!roiItem) continue;
            
            Sint32 refROINumber = 0;
            roiItem->findAndGetSint32(DCM_ReferencedROINumber, refROINumber);
            
            auto structure = std::make_unique<Structure>();
            structure->setROINumber(refROINumber);
            structure->setName(roiNames.count(refROINumber) ? roiNames[refROINumber] : "Unknown");
            
            // Get color
            const Uint16* colorData = nullptr;
            unsigned long colorCount = 0;
            if (roiItem->findAndGetUint16Array(DCM_ROIDisplayColor, colorData, &colorCount).good() && colorCount >= 3) {
                structure->setColor(
                    static_cast<uint8_t>(colorData[0]),
                    static_cast<uint8_t>(colorData[1]),
                    static_cast<uint8_t>(colorData[2])
                );
            }
            
            // Determine type from name
            std::string name = structure->getName();
            std::string upperName = name;
            std::transform(upperName.begin(), upperName.end(), upperName.begin(), ::toupper);
            
            if (upperName.find("PTV") != std::string::npos || 
                upperName.find("GTV") != std::string::npos ||
                upperName.find("CTV") != std::string::npos) {
                structure->setType("TARGET");
                structure->setPriority(1);
            } else if (upperName.find("BODY") != std::string::npos ||
                       upperName.find("EXTERNAL") != std::string::npos) {
                structure->setType("EXTERNAL");
                structure->setPriority(5);
            } else {
                structure->setType("OAR");
                structure->setPriority(3);
            }
            
            // Read contours
            DcmSequenceOfItems* contourSeq = nullptr;
            if (roiItem->findAndGetSequence(DCM_ContourSequence, contourSeq).good()) {
                for (unsigned long c = 0; c < contourSeq->card(); ++c) {
                    DcmItem* contourItem = contourSeq->getItem(c);
                    if (!contourItem) continue;
                    
                    Contour contour;
                    
                    Sint32 numPoints = 0;
                    contourItem->findAndGetSint32(DCM_NumberOfContourPoints, numPoints);
                    
                    const Float64* contourData = nullptr;
                    unsigned long dataCount = 0;
                    if (contourItem->findAndGetFloat64Array(DCM_ContourData, contourData, &dataCount).good() && contourData) {
                        size_t numCoords = dataCount / 3;
                        contour.points.reserve(numCoords);
                        
                        for (size_t p = 0; p < numCoords; ++p) {
                            contour.points.push_back({
                                contourData[p * 3],
                                contourData[p * 3 + 1],
                                contourData[p * 3 + 2]
                            });
                        }
                        
                        if (!contour.points.empty()) {
                            contour.zPosition = contour.points[0][2];
                        }
                    }
                    
                    if (!contour.points.empty()) {
                        structure->addContour(contour);
                    }
                }
            }
            
            Logger::info("  - " + structure->getName() + " (" + structure->getType() + 
                        ", " + std::to_string(structure->getContourCount()) + " contours)");
            structures->addStructure(std::move(structure));
        }
    }
    
    return structures;
#else
    Logger::warn("DCMTK not available");
    return nullptr;
#endif
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
