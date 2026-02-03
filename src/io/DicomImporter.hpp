#pragma once

#include "IDataImporter.hpp"
#include <vector>
#include <filesystem>

namespace optirad {

class DicomImporter : public IDataImporter {
public:
    bool canImport(const std::string& path) const override;
    std::unique_ptr<Patient> importPatient(const std::string& path) override;
    std::unique_ptr<StructureSet> importStructures(const std::string& path) override;
    
    // DICOM-specific methods
    bool loadDirectory(const std::string& dirPath);
    bool loadCTSeries(const std::string& dirPath);
    bool loadRTStruct(const std::string& filePath);
    bool loadRTPlan(const std::string& filePath);
    bool loadRTDose(const std::string& filePath);
    
private:
    std::vector<std::filesystem::path> m_ctFiles;
    std::filesystem::path m_rtStructFile;
    std::filesystem::path m_rtPlanFile;
    std::filesystem::path m_rtDoseFile;
    
    void scanDirectory(const std::filesystem::path& dir);
    std::string getSOPClassUID(const std::filesystem::path& file) const;
};

} // namespace optirad
