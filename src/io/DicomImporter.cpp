#include "DicomImporter.hpp"
#include "core/Patient.hpp"
#include "geometry/StructureSet.hpp"
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
        // Check if directory contains any files (DICOM often has no extension)
        int fileCount = 0;
        for (const auto& entry : fs::recursive_directory_iterator(path)) {
            if (entry.is_regular_file()) {
                auto filename = entry.path().filename().string();
                // Skip macOS metadata files
                if (filename.substr(0, 2) == "._") continue;
                fileCount++;
            }
        }
        Logger::info("Found " + std::to_string(fileCount) + " files in directory");
        return fileCount > 0;
    }
    
    // Single file
    return fs::is_regular_file(path);
}

std::unique_ptr<Patient> DicomImporter::importPatient(const std::string& path) {
    if (!loadDirectory(path)) {
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
    auto structures = std::make_unique<StructureSet>();
    
#ifdef OPTIRAD_HAS_DCMTK
    if (m_rtStructFile.empty()) {
        loadDirectory(path);
    }
    
    if (!m_rtStructFile.empty()) {
        Logger::info("Loading RT Structure Set: " + m_rtStructFile.string());
        
        DcmFileFormat fileFormat;
        OFCondition status = fileFormat.loadFile(m_rtStructFile.c_str());
        
        if (status.good()) {
            DcmDataset* dataset = fileFormat.getDataset();
            
            // Find the StructureSetROISequence
            DcmItem* roiSeqItem = nullptr;
            DcmSequenceOfItems* roiSequence = nullptr;
            
            if (dataset->findAndGetSequence(DCM_StructureSetROISequence, roiSequence).good()) {
                unsigned long numROIs = roiSequence->card();
                Logger::info("Found " + std::to_string(numROIs) + " ROIs in structure set");
                
                for (unsigned long i = 0; i < numROIs; ++i) {
                    roiSeqItem = roiSequence->getItem(i);
                    if (!roiSeqItem) continue;
                    
                    OFString roiName, roiNumber;
                    roiSeqItem->findAndGetOFString(DCM_ROIName, roiName);
                    roiSeqItem->findAndGetOFString(DCM_ROINumber, roiNumber);
                    
                    // Create structure
                    auto structure = std::make_unique<Structure>();
                    structure->setName(roiName.c_str());
                    
                    // Try to determine type from name
                    std::string name = roiName.c_str();
                    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
                    
                    if (name.find("PTV") != std::string::npos || 
                        name.find("GTV") != std::string::npos ||
                        name.find("CTV") != std::string::npos) {
                        structure->setType("TARGET");
                    } else if (name.find("LUNG") != std::string::npos ||
                               name.find("HEART") != std::string::npos ||
                               name.find("SPINAL") != std::string::npos ||
                               name.find("ESOPHAGUS") != std::string::npos) {
                        structure->setType("OAR");
                    } else if (name.find("BODY") != std::string::npos ||
                               name.find("EXTERNAL") != std::string::npos) {
                        structure->setType("EXTERNAL");
                    } else {
                        structure->setType("UNKNOWN");
                    }
                    
                    Logger::info("  - " + std::string(roiName.c_str()) + " (" + structure->getType() + ")");
                    structures->addStructure(std::move(structure));
                }
            } else {
                Logger::warn("No StructureSetROISequence found in RT Structure Set");
            }
        } else {
            Logger::error("Failed to load RT Structure Set file: " + std::string(status.text()));
        }
    }
#else
    Logger::warn("DCMTK not available - returning empty structure set");
#endif
    
    return structures;
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
    
    // Sort CT files by path (instance numbers are in path for this dataset)
    std::sort(m_ctFiles.begin(), m_ctFiles.end());
    
    Logger::info("Scan complete:");
    Logger::info("  CT slices: " + std::to_string(m_ctFiles.size()));
    if (!m_rtStructFile.empty()) Logger::info("  RT Structure Set: " + m_rtStructFile.filename().string());
    if (!m_rtPlanFile.empty()) Logger::info("  RT Plan: " + m_rtPlanFile.filename().string());
    if (!m_rtDoseFile.empty()) Logger::info("  RT Dose: " + m_rtDoseFile.filename().string());
    
    return !m_ctFiles.empty();
}

void DicomImporter::scanDirectory(const std::filesystem::path& dir) {
    namespace fs = std::filesystem;
    
    int filesScanned = 0;
    
    for (const auto& entry : fs::recursive_directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        
        // Skip hidden files and macOS metadata
        auto filename = entry.path().filename().string();
        if (filename[0] == '.') continue;
        if (filename.substr(0, 2) == "._") continue;
        if (filename == "DICOMDIR") continue;
        
        auto ext = entry.path().extension().string();
        std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
        
        // Skip obvious non-DICOM extensions
        if (ext == ".txt" || ext == ".log" || ext == ".xml" || ext == ".json") {
            continue;
        }
        
        filesScanned++;
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
    
    Logger::debug("Scanned " + std::to_string(filesScanned) + " files");
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
