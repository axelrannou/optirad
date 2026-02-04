#pragma once

#include "IPanel.hpp"
#include "core/PatientData.hpp"
#include "io/DicomImporter.hpp"
#include <memory>
#include <string>

namespace optirad {

class PatientPanel : public IPanel {
public:
    PatientPanel();
    
    void render() override;
    void update() override;  // Declare but don't define inline
    std::string getName() const override { return "Patient"; }
    
    // Access to loaded data
    PatientData* getPatientData() { return m_patientData.get(); }
    bool hasData() const { return m_patientData != nullptr; }
    
private:
    void renderImportDialog();
    void renderPatientInfo();
    void renderStructureList();
    
    void importDicom(const std::string& path);
    
    std::unique_ptr<PatientData> m_patientData;
    DicomImporter m_importer;
    
    // UI state
    bool m_showImportDialog = false;
    char m_dicomPath[512] = "";
    bool m_isImporting = false;
};

} // namespace optirad
