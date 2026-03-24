#pragma once

#include "IDataImporter.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include <vector>
#include <filesystem>
#include <memory>
#include <utility>

namespace optirad {

class PatientData;
template<typename T> class Volume;

class DicomImporter : public IDataImporter {
public:
    bool canImport(const std::string& path) const override;
    std::unique_ptr<Patient> importPatient(const std::string& path) override;
    std::unique_ptr<StructureSet> importStructures(const std::string& path) override;
    
    // Full import - returns complete PatientData
    std::unique_ptr<PatientData> importAll(const std::string& dirPath);
    
    // DICOM-specific methods
    bool loadDirectory(const std::string& dirPath);
    bool loadCTSeries(const std::string& dirPath);
    bool loadRTStruct(const std::string& filePath);
    bool loadRTPlan(const std::string& filePath);
    bool loadRTDose(const std::string& filePath);
    
    // Import CT volume with actual voxel data
    std::unique_ptr<Volume<int16_t>> importCTVolume();
    
    // Import structures with contour geometry
    std::unique_ptr<StructureSet> importStructuresWithContours();

    // Import RT Dose (returns dose matrix + its grid, or nullptr if not available)
    std::pair<std::shared_ptr<DoseMatrix>, std::shared_ptr<Grid>> importRTDose();
    
private:
    std::vector<std::filesystem::path> m_ctFiles;
    std::filesystem::path m_rtStructFile;
    std::filesystem::path m_rtPlanFile;
    std::filesystem::path m_rtDoseFile;
    
    void scanDirectory(const std::filesystem::path& dir);
    std::string getSOPClassUID(const std::filesystem::path& file) const;
    void sortCTFilesByPosition();
};

} // namespace optirad
