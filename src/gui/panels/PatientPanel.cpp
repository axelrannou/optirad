#include "PatientPanel.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <filesystem>

namespace optirad {

PatientPanel::PatientPanel() {
    // Initialize with user's home directory or last used path
    const char* home = std::getenv("HOME");
    if (home) {
        strncpy(m_dicomPath, home, sizeof(m_dicomPath) - 1);
    }
}

void PatientPanel::render() {
    ImGui::Begin("Patient Data");
    
    // Import button
    if (ImGui::Button("Import DICOM Directory", ImVec2(-1, 0))) {
        m_showImportDialog = true;
    }
    
    ImGui::Separator();
    
    // Show patient info if loaded
    if (m_patientData) {
        renderPatientInfo();
        ImGui::Separator();
        renderStructureList();
    } else {
        ImGui::TextDisabled("No patient data loaded");
    }
    
    // Import dialog
    if (m_showImportDialog) {
        renderImportDialog();
    }
    
    ImGui::End();
}

void PatientPanel::renderImportDialog() {
    ImGui::OpenPopup("Import DICOM");
    
    ImVec2 center = ImGui::GetMainViewport()->GetCenter();
    ImGui::SetNextWindowPos(center, ImGuiCond_Appearing, ImVec2(0.5f, 0.5f));
    
    if (ImGui::BeginPopupModal("Import DICOM", &m_showImportDialog, ImGuiWindowFlags_AlwaysAutoResize)) {
        ImGui::Text("Select DICOM directory:");
        ImGui::Spacing();
        
        ImGui::InputText("Path", m_dicomPath, sizeof(m_dicomPath));
        ImGui::SameLine();
        if (ImGui::Button("Browse...")) {
            // TODO: Add file browser (or use native dialog)
            Logger::info("File browser not yet implemented");
        }
        
        ImGui::Spacing();
        
        if (m_isImporting) {
            ImGui::TextColored(ImVec4(1.0f, 1.0f, 0.0f, 1.0f), "Importing...");
            ImGui::ProgressBar(-1.0f * ImGui::GetTime());
        } else {
            if (ImGui::Button("Import", ImVec2(120, 0))) {
                importDicom(m_dicomPath);
            }
            ImGui::SameLine();
            if (ImGui::Button("Cancel", ImVec2(120, 0))) {
                m_showImportDialog = false;
            }
        }
        
        ImGui::EndPopup();
    }
}

void PatientPanel::renderPatientInfo() {
    if (!m_patientData) return;
    
    ImGui::Text("Patient Information");
    ImGui::Indent();
    
    if (auto* patient = m_patientData->getPatient()) {
        ImGui::Text("Name: %s", patient->getName().c_str());
        ImGui::Text("ID:   %s", patient->getID().c_str());
    }
    
    if (auto* ct = m_patientData->getCTVolume()) {
        auto dims = ct->getGrid().getDimensions();
        auto spacing = ct->getGrid().getSpacing();
        
        ImGui::Text("CT Volume:");
        ImGui::Text("  Dimensions: %zux%zux%zu", dims[0], dims[1], dims[2]);
        ImGui::Text("  Spacing: %.2fx%.2fx%.2f mm", spacing[0], spacing[1], spacing[2]);
        ImGui::Text("  Voxels: %zu", ct->size());
    }
    
    ImGui::Unindent();
}

void PatientPanel::renderStructureList() {
    if (!m_patientData) return;
    
    auto* structures = m_patientData->getStructureSet();
    if (!structures || structures->getCount() == 0) return;
    
    ImGui::Text("Structures (%zu)", structures->getCount());
    
    if (ImGui::BeginTable("Structures", 4, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg)) {
        ImGui::TableSetupColumn("Visible");
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("Type");
        ImGui::TableSetupColumn("Contours");
        ImGui::TableHeadersRow();
        
        for (size_t i = 0; i < structures->getCount(); ++i) {
            auto* s = structures->getStructure(i);
            if (!s) continue;
            
            ImGui::TableNextRow();
            
            // Visible checkbox
            ImGui::TableNextColumn();
            bool visible = s->isVisible();
            if (ImGui::Checkbox(("##vis" + std::to_string(i)).c_str(), &visible)) {
                const_cast<Structure*>(s)->setVisible(visible);
            }
            
            // Name with color indicator
            ImGui::TableNextColumn();
            auto color = s->getColor();
            ImGui::ColorButton(("##color" + std::to_string(i)).c_str(),
                             ImVec4(color[0]/255.0f, color[1]/255.0f, color[2]/255.0f, 1.0f),
                             ImGuiColorEditFlags_NoTooltip | ImGuiColorEditFlags_NoBorder);
            ImGui::SameLine();
            ImGui::Text("%s", s->getName().c_str());
            
            // Type
            ImGui::TableNextColumn();
            ImGui::Text("%s", s->getType().c_str());
            
            // Contour count
            ImGui::TableNextColumn();
            ImGui::Text("%zu", s->getContourCount());
        }
        
        ImGui::EndTable();
    }
}

void PatientPanel::importDicom(const std::string& path) {
    namespace fs = std::filesystem;
    
    if (!fs::exists(path)) {
        Logger::error("Path does not exist: " + path);
        return;
    }
    
    m_isImporting = true;
    Logger::info("Importing DICOM from: " + path);
    
    // Import in main thread (for now - TODO: move to background thread)
    m_patientData = m_importer.importAll(path);
    
    if (m_patientData) {
        Logger::info("DICOM import successful");
        m_showImportDialog = false;
    } else {
        Logger::error("DICOM import failed");
    }
    
    m_isImporting = false;
}

void PatientPanel::update() {
    // No per-frame update needed
}

} // namespace optirad
