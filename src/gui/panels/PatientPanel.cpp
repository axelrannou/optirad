#include "PatientPanel.hpp"
#include "../Layout.hpp"
#include "io/DicomExporter.hpp"
#include "utils/Logger.hpp"
#include "DoseStatsCsv.hpp"
#include "ImportedPlanDoseView.hpp"
#include <imgui.h>
#include <filesystem>
#include <algorithm>
#include <fstream>
#include <string>
#include <unordered_map>
#include <chrono>
#include <stdexcept>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad {

PatientPanel::PatientPanel(GuiAppState& state, unsigned int eyeOpenTexture,
                           unsigned int eyeClosedTexture)
    : m_state(state), m_eyeOpenTexture(eyeOpenTexture), m_eyeClosedTexture(eyeClosedTexture) {
    // Default to project data/ directory
    snprintf(m_dicomPath, sizeof(m_dicomPath), "%s", OPTIRAD_DATA_DIR);
    m_dicomPath[sizeof(m_dicomPath) - 1] = '\0';
    // Default export path
    snprintf(m_exportPath, sizeof(m_exportPath), "%s/export", OPTIRAD_DATA_DIR);
    m_exportPath[sizeof(m_exportPath) - 1] = '\0';
    m_state.exportDir = m_exportPath;
}

PatientPanel::~PatientPanel() {
    if (m_importFuture.valid()) m_importFuture.wait();
    if (m_isImporting) m_state.taskRunning = false;
}

void PatientPanel::render() {
    update();
    if (!m_visible) return;

    m_state.exportDir = m_exportPath;

    ImGui::Begin("Patient Data", &m_visible, ImGuiWindowFlags_HorizontalScrollbar);

    // Equal-width buttons stack when the labels cannot fit side by side.
    {
        const float spacing = ImGui::GetStyle().ItemSpacing.x;
        const DoseEntry* selEntry = m_state.doseStore.getSelected();
        const bool hasDose = selEntry && selEntry->dose;
        const bool hasSeq  = hasDose && m_state.seqCache.count(selEntry->id) > 0;
        const char* exportLabel = hasSeq ? "Export Dose + RT Plan" : "Export Dose";
        const float minWidth = std::max(layout::framedTextWidth("Import DICOM"),
                                        layout::framedTextWidth(exportLabel));
        const bool sideBySide = ImGui::GetContentRegionAvail().x >= minWidth * 2.0f + spacing;
        const float btnW = std::max(minWidth, sideBySide ?
            (ImGui::GetContentRegionAvail().x - spacing) * 0.5f : ImGui::GetContentRegionAvail().x);

        ImGui::BeginDisabled(m_state.taskRunning);
        if (ImGui::Button("Import DICOM", ImVec2(btnW, 0)))
            m_showImportDialog = true;
        ImGui::EndDisabled();

        if (sideBySide) ImGui::SameLine();

        ImGui::BeginDisabled(!hasDose || m_isImporting);
        if (ImGui::Button(exportLabel, ImVec2(btnW, 0))) {
            const std::string outDir =
                (m_exportPath[0] != '\0') ? m_exportPath : (OPTIRAD_DATA_DIR "/export");
            DicomExporter exporter;

            // Export the RT Plan first so its SOP Instance UID can be cross-referenced
            // from the RT Dose (required by the RT Dose IOD when DoseSummationType is PLAN).
            std::string planSOPInstanceUID;
            bool ok = true;
            if (hasSeq) {
                const auto& seqEntry = m_state.seqCache.at(selEntry->id);
                ok = exporter.exportRTPlan(*m_state.plan, *m_state.stf,
                                           seqEntry.sequences, m_state.dicomContext, outDir,
                                           &planSOPInstanceUID);
            }
            if (ok) {
                exporter.exportRTDose(*selEntry->dose, m_state.dicomContext,
                                       outDir, selEntry->name, planSOPInstanceUID);
            }
        }
        ImGui::EndDisabled();
    }

    ImGui::Separator();
    renderImportProgress();
    
    // Show patient info if loaded
    if (m_state.patientData) {
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
    
    if (ImGui::BeginPopupModal("Import DICOM", &m_showImportDialog, ImGuiWindowFlags_HorizontalScrollbar)) {
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
        
        ImGui::BeginDisabled(m_state.taskRunning);
        if (ImGui::Button("Import", ImVec2(120, 0))) {
            importDicom(m_dicomPath);
            if (m_isImporting) ImGui::CloseCurrentPopup();
        }
        ImGui::EndDisabled();
        layout::sameLineIfFits(120.0f);
        if (ImGui::Button("Cancel", ImVec2(120, 0))) {
            m_showImportDialog = false;
        }

        ImGui::EndPopup();
    }
}

void PatientPanel::renderImportProgress() {
    if (m_isImporting) {
        ImGui::TextWrapped("%s", m_state.taskStatus.c_str());
        const char* stages[] = {"CT + dose", "Structures", "RTPLAN"};
        const int active = static_cast<int>(m_importStage);
        for (int i = 0; i < 3; ++i) {
            ImGui::Text("%s", stages[i]);
            const float progress = i < active ? 1.0f :
                (i == active ? -static_cast<float>(ImGui::GetTime()) : 0.0f);
            ImGui::ProgressBar(progress, ImVec2(-1, 0),
                               i < active ? "Loaded" : (i == active ? "Loading..." : "Waiting"));
        }
        ImGui::Separator();
    }
    if (!m_importError.empty()) {
        ImGui::TextColored(ImVec4(1.0f, 0.4f, 0.4f, 1.0f), "Import error:");
        ImGui::TextWrapped("%s", m_importError.c_str());
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
    ImGui::BeginChild("FileBrowser", ImVec2(-1, browserH), true, ImGuiWindowFlags_HorizontalScrollbar);
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
            if (ImGui::Selectable(label.c_str(), false, 0,
                                  ImVec2(ImGui::CalcTextSize(label.c_str()).x, 0))) {
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
    if (!m_state.patientData) return;
    
    if (ImGui::CollapsingHeader("Patient Information", ImGuiTreeNodeFlags_DefaultOpen)) {
        ImGui::Indent();
        
        if (auto* patient = m_state.patientData->getPatient()) {
            ImGui::Text("Name: %s", patient->getName().c_str());
            ImGui::Text("ID:   %s", patient->getID().c_str());
        }
        
        ImGui::Unindent();
    }
    
    if (auto* ct = m_state.patientData->getCTVolume()) {
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
    if (!m_state.patientData) return;
    
    auto* structures = m_state.patientData->getStructureSet();
    if (!structures || structures->getCount() == 0) return;
    
    bool allVisible = true;
    for (size_t i = 0; i < structures->getCount(); ++i) {
        const auto* structure = structures->getStructure(i);
        if (structure && !structure->isVisible()) {
            allVisible = false;
            break;
        }
    }

    auto visibilityButton = [&](const char* id, bool visible, const char* tooltip) {
        const unsigned int texture = visible ? m_eyeOpenTexture : m_eyeClosedTexture;
        ImGui::PushStyleColor(ImGuiCol_Button, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive, ImVec4(0, 0, 0, 0));
        ImGui::PushStyleVar(ImGuiStyleVar_FrameBorderSize, 0.0f);
        bool clicked;
        if (texture) {
            const float size = ImGui::GetTextLineHeight();
            clicked = ImGui::ImageButton(id,
                reinterpret_cast<void*>(static_cast<intptr_t>(texture)),
                ImVec2(size * 1.58f, size), ImVec2(0, 0), ImVec2(1, 1),
                ImVec4(0, 0, 0, 0), ImGui::GetStyleColorVec4(ImGuiCol_Text));
        } else {
            clicked = ImGui::SmallButton((std::string(visible ? "Hide" : "Show") + id).c_str());
        }
        ImGui::PopStyleVar();
        ImGui::PopStyleColor(3);
        if (ImGui::IsItemHovered())
            ImGui::SetTooltip("%s %s", visible ? "Hide" : "Show", tooltip);
        return clicked;
    };

    if (ImGui::BeginTable("Structures", 3, layout::contentTableFlags,
                          layout::tableSize(structures->getCount(), ImGui::GetFrameHeight(),
                                            ImGui::GetFrameHeight()))) {
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("Type");
        ImGui::TableSetupColumn("Contours");

        ImGui::TableNextRow(ImGuiTableRowFlags_Headers);
        ImGui::TableNextColumn();
        if (visibilityButton("##allVisibility", allVisible, "all structures")) {
            for (size_t i = 0; i < structures->getCount(); ++i) {
                if (auto* structure = structures->getStructure(i))
                    structure->setVisible(!allVisible);
            }
        }
        ImGui::SameLine();
        ImGui::TableHeader("Name");
        ImGui::TableNextColumn();
        ImGui::TableHeader("Type");
        ImGui::TableNextColumn();
        ImGui::TableHeader("Contours");
        
        for (size_t i = 0; i < structures->getCount(); ++i) {
            auto* s = structures->getStructure(i);
            if (!s) continue;
            
            ImGui::TableNextRow();
            
            ImGui::TableNextColumn();
            const bool visible = s->isVisible();
            const std::string id = "##vis" + std::to_string(i);
            if (visibilityButton(id.c_str(), visible, s->getName().c_str()))
                s->setVisible(!visible);
            
            // Name with color indicator
            ImGui::SameLine();
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

    ImGui::BeginDisabled(m_isImporting);
    if (ImGui::Button("Export All Stats (CSV)")) {
        exportAllStatsCSV();
    }
    ImGui::EndDisabled();
    ImGui::Spacing();

    if (ImGui::BeginTable("DoseList", 5,
            ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg | layout::contentTableFlags,
            layout::tableSize(static_cast<size_t>(dm.count()), ImGui::GetFrameHeight()))) {
        ImGui::TableSetupColumn("##sel", ImGuiTableColumnFlags_WidthFixed, ImGui::GetFrameHeight());
        ImGui::TableSetupColumn("Name");
        ImGui::TableSetupColumn("Performed on");
        ImGui::TableSetupColumn("Max (Gy)");
        ImGui::TableSetupColumn("##del");
        ImGui::TableHeadersRow();

        bool taskBusy = m_state.taskRunning;

        std::unordered_map<int, const DoseEntry*> entriesById;
        entriesById.reserve(static_cast<size_t>(dm.count()));
        for (const auto& entry : dm.getEntries()) {
            entriesById.emplace(entry.id, &entry);
        }

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
            float colMaxX = colMinX + ImGui::GetContentRegionAvail().x;

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
                ImGui::SetNextItemWidth(layout::framedTextWidth(m_editDoseName));
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
                    auto linkedIt = entriesById.find(linkedId);
                    if (linkedIt != entriesById.end()) {
                        ImGui::Text("%s", linkedIt->second->name.c_str());
                    } else {
                        ImGui::TextDisabled("-");
                    }
                } else {
                    ImGui::TextDisabled("-");
                }
            }

            // Max dose
            ImGui::TableNextColumn();
            if (entry->dose) {
                ImGui::Text("%.2f", entry->maxDoseGy);
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

    if (m_state.taskRunning) {
        Logger::error("Cannot import DICOM while a task is running");
        return;
    }
    
    std::error_code ec;
    if (!fs::is_directory(path, ec)) {
        m_importError = "Invalid DICOM directory: " + path;
        if (ec) m_importError += " (" + ec.message() + ")";
        Logger::error(m_importError);
        return;
    }
    
    m_importPath = path;
    m_importError.clear();
    m_isImporting = true;
    m_state.taskRunning = true;
    m_importStage = ImportStage::Images;
    m_importedDoseId = -1;
    m_showImportDialog = false;
    Logger::info("Importing DICOM from: " + path);
    try {
        startImportStage();
    } catch (const std::exception& e) {
        m_importError = "DICOM import failed: " + std::string(e.what());
        Logger::error(m_importError);
        m_isImporting = false;
        m_state.taskRunning = false;
        m_state.taskStatus = m_importError;
    }
}

void PatientPanel::setImportStatus(const std::string& status) {
    std::lock_guard<std::mutex> lock(m_importMutex);
    m_importStatus = status;
}

void PatientPanel::startImportStage() {
    const auto stage = m_importStage;
    setImportStatus(stage == ImportStage::Images ? "Scanning DICOM directory" :
                    (stage == ImportStage::Structures ? "Parsing RTSTRUCT contours" : "Loading RTPLAN"));
    const auto patient = m_state.patientData;
    m_importFuture = std::async(std::launch::async, [this, stage, patient]() {
        ImportResult result;
        auto status = [this](const std::string& text) { setImportStatus(text); };
        if (stage == ImportStage::Images) {
            result.patient = m_importer.importImages(m_importPath, status);
            if (!result.patient || !result.patient->hasValidCT())
                throw std::runtime_error("Failed to load a valid CT from " + m_importPath);
        } else if (stage == ImportStage::Structures) {
            result.structures = m_importer.importPreparedStructures(patient->getCTVolume(), status);
        } else if (stage == ImportStage::Plan) {
            result.plan = m_importer.importRTPlan();
            if (result.plan) {
                for (const auto& warning : result.plan->warnings) Logger::warn("RTPLAN: " + warning);
                status("Building RTPLAN beams and delivery");
                try {
                    result.preparedPlan = ImportedPlanDoseView::prepareImported(patient, *result.plan);
                } catch (const std::exception& e) {
                    result.planError = "RTPLAN plan could not be built: " + std::string(e.what());
                    Logger::error(result.planError);
                }
            }
        }
        return result;
    });
    m_state.taskProgress = static_cast<float>(stage) / 3.0f;
    std::lock_guard<std::mutex> lock(m_importMutex);
    m_state.taskStatus = m_importStatus;
}

void PatientPanel::exportAllStatsCSV() {
    namespace fs = std::filesystem;

    if (!m_state.patientData) {
        Logger::error("DoseStats CSV export: no patient data loaded");
        return;
    }

    auto& dm = m_state.doseStore;
    if (dm.count() == 0) return;

    std::string patientID;
    if (auto* patient = m_state.patientData->getPatient()) patientID = patient->getID();
    if (patientID.empty()) patientID = "UnknownPatient";

    const std::string outDir =
        (!m_state.exportDir.empty()) ? m_state.exportDir : (OPTIRAD_DATA_DIR "/export");
    if (!fs::exists(outDir)) {
        std::error_code ec;
        if (!fs::create_directories(outDir, ec)) {
            Logger::error("DoseStats CSV export: cannot create output directory: " + outDir
                          + " (" + ec.message() + ")");
            return;
        }
    }

    const std::string filename =
        (fs::path(outDir) / ("DoseStats_" + sanitizeFilename(patientID) + "_All.csv")).string();
    std::ofstream f(filename);
    if (!f) {
        Logger::error("DoseStats CSV export: cannot open file for writing: " + filename);
        return;
    }

    double rxDose = (m_state.plan && m_state.plan->getPrescribedDose() > 0)
                    ? m_state.plan->getPrescribedDose() : 60.0;

    writeDoseStatsCsvHeader(f);
    for (int i = 0; i < dm.count(); ++i) {
        const auto* entry = dm.getEntry(i);
        if (!entry || !entry->dose || !entry->grid) continue;
        const auto& stats = m_state.getOrComputeStats(i, *m_state.patientData, rxDose);
        writeDoseStatsCsvRows(f, entry->name, stats);
    }

    Logger::info("DoseStats CSV export: written to " + filename);
}

void PatientPanel::update() {
    if (!m_isImporting) return;
    {
        std::lock_guard<std::mutex> lock(m_importMutex);
        m_state.taskStatus = m_importStatus;
    }
    if (!m_importFuture.valid() ||
        m_importFuture.wait_for(std::chrono::seconds(0)) != std::future_status::ready) return;

    try {
        auto result = m_importFuture.get();
        if (m_importStage == ImportStage::Images) {
            m_state.resetAllDoses();
            m_state.resetPlan();
            m_state.resetLeafSequence();
            m_state.patientData = std::move(result.patient);
            m_editingDoseIdx = -1;
            m_state.dicomContext = m_importer.getContext();
            auto [dose, grid] = m_importer.takeImportedDose();
            if (dose) m_importedDoseId = m_state.doseStore.addDose("Imported (DICOM)", dose, grid);
            m_state.syncSelectedDose();
            m_importStage = ImportStage::Structures;
        } else if (m_importStage == ImportStage::Structures) {
            m_state.patientData->setStructureSet(std::move(result.structures));
            m_state.statsCache.clear();
            m_state.dicomContext = m_importer.getContext();
            m_importStage = ImportStage::Plan;
        } else {
            m_state.patientData->setImportedPlan(std::move(result.plan));
            if (result.preparedPlan.built.plan)
                ImportedPlanDoseView::installPrepared(m_state, m_importedDoseId, std::move(result.preparedPlan));
            m_importError = std::move(result.planError);
            m_importStage = ImportStage::Complete;
            m_isImporting = false;
            m_state.taskRunning = false;
            m_state.taskProgress = 1.0f;
            m_state.taskStatus = m_importError.empty() ? "DICOM import complete" : m_importError;
            Logger::info(m_importError.empty() ? "DICOM import successful" : "DICOM images and structures imported; RTPLAN build failed");
        }
        if (m_isImporting) startImportStage();
    } catch (const std::exception& e) {
        m_importError = "DICOM import failed: " + std::string(e.what());
        Logger::error(m_importError);
        m_isImporting = false;
        m_state.taskRunning = false;
        m_state.taskStatus = m_importError;
    }
}

} // namespace optirad
