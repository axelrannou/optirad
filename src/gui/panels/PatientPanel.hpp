#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"
#include "core/PatientData.hpp"
#include "io/DicomImporter.hpp"
#include "ImportedPlanDoseView.hpp"
#include <future>
#include <memory>
#include <mutex>
#include <string>

namespace optirad {

class PatientPanel : public IPanel {
public:
    PatientPanel(GuiAppState& state, unsigned int eyeOpenTexture = 0,
                 unsigned int eyeClosedTexture = 0);
    ~PatientPanel() override;
    
    void render() override;
    void update() override;
    std::string getName() const override { return "Patient"; }
    
    // Access to loaded data
    PatientData* getPatientData() { return m_state.patientData.get(); }
    bool hasData() const { return m_state.patientData != nullptr; }
    bool isImporting() const { return m_isImporting; }

    void importDicom(const std::string& path);
    
private:
    void renderImportDialog();
    void renderFileBrowser();
    void renderPatientInfo();
    void renderStructureList();
    void renderDoseList();
    void renderImportProgress();
    void startImportStage();
    void setImportStatus(const std::string& status);
    
    void exportAllStatsCSV();
    
    GuiAppState& m_state;
    DicomImporter m_importer;
    unsigned int m_eyeOpenTexture;
    unsigned int m_eyeClosedTexture;
    
    // UI state
    bool m_showImportDialog = false;
    char m_dicomPath[512] = "";
    bool m_isImporting = false;
    enum class ImportStage { Images, Structures, Plan, Complete };
    ImportStage m_importStage = ImportStage::Complete;
    struct ImportResult {
        std::unique_ptr<PatientData> patient;
        std::unique_ptr<StructureSet> structures;
        std::shared_ptr<ImportedPlan> plan;
        ImportedPlanDoseView::PreparedPlan preparedPlan;
        std::string planError;
    };
    std::future<ImportResult> m_importFuture;
    std::mutex m_importMutex;
    std::string m_importStatus;
    std::string m_importError;
    std::string m_importPath;
    int m_importedDoseId = -1;

    // Dose name inline editing
    int m_editingDoseIdx = -1;
    char m_editDoseName[128] = "";

    // RT Plan export
    char m_exportPath[512] = "";
};

} // namespace optirad
