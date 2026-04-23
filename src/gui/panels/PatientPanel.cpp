#include "PatientPanel.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <filesystem>
#include <algorithm>
#include <string>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad {

PatientPanel::PatientPanel(GuiAppState& state) : m_state(state) {
    // Default to project data/ directory
    snprintf(m_dicomPath, sizeof(m_dicomPath), "%s", OPTIRAD_DATA_DIR);
    m_dicomPath[sizeof(m_dicomPath) - 1] = '\0';
}

void PatientPanel::render() {
    if (!m_visible) return;

    ImGui::Begin("Patient Data", &m_visible);
    
    // Import button
    if (ImGui::Button("Import DICOM Directory", ImVec2(-1, 0))) {
        m_showImportDialog = true;
    }
    
    ImGui::Separator();
    
    // Show patient info if loaded
    if (m_patientData) {
        renderPatientInfo();
        ImGui::Separator();

        if (ImGui::CollapsingHeader("Dose", ImGuiTreeNodeFlags_DefaultOpen)) {
            renderDoseList();
        }
        
        if (ImGui::CollapsingHeader("Structures", ImGuiTreeNodeFlags_DefaultOpen)) {
            renderStructureList();
        }
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
    ImGui::SetNextWindowSize(ImVec2(600, 500), ImGuiCond_Appearing);
    
    if (ImGui::BeginPopupModal("Import DICOM", &m_showImportDialog)) {
        ImGui::Text("Select DICOM directory:");
        ImGui::Spacing();
        
        ImGui::InputText("Path", m_dicomPath, sizeof(m_dicomPath));
        
        ImGui::Spacing();
        ImGui::Separator();
        ImGui::Spacing();
        
        // Inline file browser
        renderFileBrowser();
        
        ImGui::Spacing();
        ImGui::Separator();
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

void PatientPanel::renderFileBrowser() {
    namespace fs = std::filesystem;

    std::string currentPath(m_dicomPath);
    if (currentPath.empty() || !fs::exists(currentPath)) {
        currentPath = OPTIRAD_DATA_DIR;
        snprintf(m_dicomPath, sizeof(m_dicomPath), "%s", currentPath.c_str());
    }

    // If path points to a file, browse its parent
    if (fs::exists(currentPath) && !fs::is_directory(currentPath)) {
        currentPath = fs::path(currentPath).parent_path().string();
    }

    // Navigate up
    if (ImGui::Button("^ Parent Directory")) {
        fs::path parent = fs::path(currentPath).parent_path();
        if (!parent.empty()) {
            snprintf(m_dicomPath, sizeof(m_dicomPath), "%s", parent.string().c_str());
        }
    }

    // Directory listing (reserve space for buttons below: spacing + separator + spacing + button row)
    float reservedH = ImGui::GetFrameHeightWithSpacing() + ImGui::GetStyle().ItemSpacing.y * 3 + 10;
    float browserH = ImGui::GetContentRegionAvail().y - reservedH;
    if (browserH < 100.0f) browserH = 100.0f;
    ImGui::BeginChild("FileBrowser", ImVec2(-1, browserH), true);
    try {
        std::vector<fs::directory_entry> dirs, files;
        for (const auto& entry : fs::directory_iterator(currentPath)) {
            if (entry.is_directory())
                dirs.push_back(entry);
            else
                files.push_back(entry);
        }
        std::sort(dirs.begin(), dirs.end(),
            [](const fs::directory_entry& a, const fs::directory_entry& b) {
                return a.path().filename() < b.path().filename();
            });
        std::sort(files.begin(), files.end(),
            [](const fs::directory_entry& a, const fs::directory_entry& b) {
                return a.path().filename() < b.path().filename();
            });

        for (const auto& dir : dirs) {
            std::string label = "[DIR]  " + dir.path().filename().string();
            if (ImGui::Selectable(label.c_str())) {
                snprintf(m_dicomPath, sizeof(m_dicomPath), "%s",
                         dir.path().string().c_str());
            }
        }
        for (const auto& file : files) {
            ImGui::TextDisabled("       %s", file.path().filename().string().c_str());
        }
    } catch (const std::exception& e) {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "Error: %s", e.what());
    }
    ImGui::EndChild();
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
        if (ImGui::CollapsingHeader("CT Volume")) {
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

void PatientPanel::renderDoseList() {
    auto& dm = m_state.doseStore;

    if (dm.count() == 0) {
        ImGui::TextDisabled("No dose data");
        return;
    }

    if (ImGui::BeginTable("DoseList", 5,
            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | ImGuiTableFlags_SizingStretchProp)) {
        ImGui::TableSetupColumn("##sel", ImGuiTableColumnFlags_WidthFixed, 30.0f);
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("Performed on", ImGuiTableColumnFlags_WidthStretch, 1.0f);
        ImGui::TableSetupColumn("Max (Gy)", ImGuiTableColumnFlags_WidthStretch, 0.5f);
        ImGui::TableSetupColumn("##del", ImGuiTableColumnFlags_WidthFixed, 25.0f);
        ImGui::TableHeadersRow();

        bool taskBusy = m_state.taskRunning;

        int toRemove = -1;
        for (int i = 0; i < dm.count(); ++i) {
            const auto* entry = dm.getEntry(i);
            if (!entry) continue;

            ImGui::TableNextRow();

            // Selectable spanning the full row as the hit area
            ImGui::TableNextColumn();
            bool selected = (i == dm.getSelectedIdx());
            if (taskBusy) ImGui::BeginDisabled();
            float rowHeight = ImGui::GetTextLineHeight() + ImGui::GetStyle().CellPadding.y * 2.0f;
            ImGui::PushStyleColor(ImGuiCol_Header,        ImVec4(0,0,0,0));
            ImGui::PushStyleColor(ImGuiCol_HeaderHovered, ImGui::GetStyleColorVec4(ImGuiCol_HeaderHovered));
            ImGui::PushStyleColor(ImGuiCol_HeaderActive,  ImGui::GetStyleColorVec4(ImGuiCol_HeaderActive));
            if (ImGui::Selectable(("##drow" + std::to_string(i)).c_str(), selected,
                    ImGuiSelectableFlags_SpanAllColumns | ImGuiSelectableFlags_AllowOverlap |
                    ImGuiSelectableFlags_AllowDoubleClick,
                    ImVec2(0, rowHeight))) {
                if (ImGui::IsMouseDoubleClicked(0)) {
                    // Double-click: enter edit mode on the name
                    m_editingDoseIdx = i;
                    snprintf(m_editDoseName, sizeof(m_editDoseName), "%s", entry->name.c_str());
                } else if (!selected) {
                    dm.setSelected(i);
                    m_state.syncSelectedDose();
                }
            }
            ImGui::PopStyleColor(3);
            if (taskBusy) ImGui::EndDisabled();

            // Draw selection indicator circle directly via DrawList
            float radius = ImGui::GetTextLineHeight() * 0.35f;
            ImVec2 cellMin = ImGui::GetItemRectMin();
            ImVec2 cellMax = ImGui::GetItemRectMax();

            float colMinX = ImGui::GetCursorScreenPos().x;
            float colMaxX = colMinX + ImGui::GetColumnWidth();

            ImVec2 center(
                (colMinX + colMaxX) * 0.5f,
                (cellMin.y + cellMax.y) * 0.5f
            );

            ImDrawList* dl = ImGui::GetWindowDrawList();
            ImU32 borderCol = ImGui::GetColorU32(ImGuiCol_Text);
            ImU32 fillCol   = ImGui::GetColorU32(selected ? ImGuiCol_Text : ImGuiCol_FrameBg);
            dl->AddCircleFilled(center, radius, fillCol);
            dl->AddCircle(center, radius, borderCol, 0, 1.2f);

            // Name (double-click to edit)
            ImGui::TableNextColumn();
            if (m_editingDoseIdx == i) {
                ImGui::SetNextItemWidth(-1);
                ImGui::SetKeyboardFocusHere();
                bool committed = ImGui::InputText(
                    ("##dname" + std::to_string(i)).c_str(),
                    m_editDoseName, sizeof(m_editDoseName),
                    ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_AutoSelectAll);
                if (committed || ImGui::IsItemDeactivatedAfterEdit() ||
                    (!ImGui::IsItemActive() && !ImGui::IsItemFocused() && ImGui::IsItemDeactivated())) {
                    if (m_editDoseName[0] != '\0') {
                        dm.renameDose(i, m_editDoseName);
                    }
                    m_editingDoseIdx = -1;
                }
                // Cancel on Escape
                if (ImGui::IsKeyPressed(ImGuiKey_Escape)) {
                    m_editingDoseIdx = -1;
                }
            } else {
                ImGui::Text("%s", entry->name.c_str());
            }

            // Performed on (show linked optimization name for leaf sequencing doses)
            ImGui::TableNextColumn();
            {
                auto seqIt = m_state.seqCache.find(entry->id);
                if (seqIt != m_state.seqCache.end()) {
                    int linkedId = seqIt->second.linkedOptDoseId;
                    bool found = false;
                    for (int j = 0; j < dm.count(); ++j) {
                        const auto* optEntry = dm.getEntry(j);
                        if (optEntry && optEntry->id == linkedId) {
                            ImGui::Text("%s", optEntry->name.c_str());
                            found = true;
                            break;
                        }
                    }
                    if (!found) ImGui::TextDisabled("-");
                } else {
                    ImGui::TextDisabled("-");
                }
            }

            // Max dose
            ImGui::TableNextColumn();
            if (entry->dose) {
                ImGui::Text("%.2f", entry->dose->getMax());
            } else {
                ImGui::TextDisabled("-");
            }

            // Delete button
            ImGui::TableNextColumn();
            if (taskBusy) ImGui::BeginDisabled();
            if (ImGui::SmallButton(("X##ddel" + std::to_string(i)).c_str())) {
                toRemove = i;
            }
            if (taskBusy) ImGui::EndDisabled();
        }

        if (toRemove >= 0) {
            dm.removeDose(toRemove);
            m_state.syncSelectedDose();
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

        // Add imported RT Dose to DoseStore if present
        auto [doseMatrix, doseGrid] = m_importer.takeImportedDose();
        if (doseMatrix) {
            m_state.doseStore.addDose("Imported (DICOM)", doseMatrix, doseGrid);
            m_state.syncSelectedDose();
        }
    } else {
        Logger::error("DICOM import failed");
    }
    
    m_isImporting = false;
}

void PatientPanel::update() {
    // No per-frame update needed
}

} // namespace optirad
