#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/PatientData.hpp"
#include "io/DicomImporter.hpp"
#include <memory>
#include <string>

namespace optirad {

class PatientPanel : public IPanel {
public:
    PatientPanel(GuiAppState& state);
    
    void render() override;
    void update() override;  // Declare but don't define inline
    std::string getName() const override { return "Patient"; }
    
    // Access to loaded data
    PatientData* getPatientData() { return m_state.patientData.get(); }
    bool hasData() const { return m_state.patientData != nullptr; }

    void importDicom(const std::string& path);
    
private:
    void renderImportDialog();
    void renderFileBrowser();
    void renderPatientInfo();
    void renderStructureList();
    void renderDoseList();
    
    void exportAllStatsCSV();
    
    GuiAppState& m_state;
    DicomImporter m_importer;
    
    // UI state
    bool m_showImportDialog = false;
    char m_dicomPath[512] = "";
    bool m_isImporting = false;

    // Dose name inline editing
    int m_editingDoseIdx = -1;
    char m_editDoseName[128] = "";

    // RT Plan export
    char m_exportPath[512] = "";
};

} // namespace optirad
