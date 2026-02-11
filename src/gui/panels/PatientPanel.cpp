#include "PatientPanel.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <filesystem>

namespace optirad {

PatientPanel::PatientPanel() {
    // Initialize with user's home directory or last used path
    const char* home = std::getenv("HOME");
    if (home) {
        snprintf(m_dicomPath, sizeof(m_dicomPath), "%s", home);
        m_dicomPath[sizeof(m_dicomPath) - 1] = '\0';  // Ensure null termination
    } else {
        m_dicomPath[0] = '\0';  // Initialize to empty string if no HOME
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
    
    if (ImGui::CollapsingHeader("Patient Information", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        
        if (auto* patient = m_patientData->getPatient()) {
            ImGui::Text("Name: %s", patient->getName().c_str());
            ImGui::Text("ID:   %s", patient->getID().c_str());
        }
        
        ImGui::Unindent();
    }
    
    if (auto* ct = m_patientData->getCTVolume()) {
        if (ImGui::CollapsingHeader("CT Volume", ImGuiTreeNodeFlags_DefaultOpen)) {
            ImGui::Indent();
            
            const auto& grid = ct->getGrid();
            auto dims = grid.getDimensions();
            auto spacing = grid.getSpacing();
            auto origin = grid.getOrigin();
            
            ImGui::Text("Dimensions: %zu x %zu x %zu", dims[0], dims[1], dims[2]);
            ImGui::Text("Spacing: %.2f x %.2f x %.2f mm", spacing[0], spacing[1], spacing[2]);
            ImGui::Text("Origin: (%.1f, %.1f, %.1f) mm", origin[0], origin[1], origin[2]);
            ImGui::Text("Total Voxels: %zu", ct->size());
            
            ImGui::Spacing();
            ImGui::Separator();
            ImGui::Text("Geometric Information:");
            
            ImGui::BulletText("Patient Position: %s", grid.getPatientPosition().c_str());
            ImGui::BulletText("Slice Thickness: %.2f mm", grid.getSliceThickness());
            
            auto orient = grid.getImageOrientation();
            ImGui::BulletText("Image Orientation (Row): [%.3f, %.3f, %.3f]", 
                             orient[0], orient[1], orient[2]);
            ImGui::BulletText("Image Orientation (Col): [%.3f, %.3f, %.3f]", 
                             orient[3], orient[4], orient[5]);
            
            // Calculate and show slice normal (cross product)
            double nx = orient[1]*orient[5] - orient[2]*orient[4];
            double ny = orient[2]*orient[3] - orient[0]*orient[5];
            double nz = orient[0]*orient[4] - orient[1]*orient[3];
            ImGui::BulletText("Slice Normal: [%.3f, %.3f, %.3f]", nx, ny, nz);
            
            ImGui::Unindent();
        }
    }
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
