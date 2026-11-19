#include <gtest/gtest.h>
#include "gui/panels/DVHPanel.hpp"
#include "gui/panels/DoseStatsPanel.hpp"
#include "gui/panels/OptimizationPanel.hpp"
#include "gui/panels/PatientPanel.hpp"
#include <imgui.h>
#include <imgui_internal.h>

namespace optirad::tests {

class PanelLayoutTest : public testing::Test {
protected:
    GuiAppState state;

    void SetUp() override {
        ImGui::CreateContext();
        auto& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(1600, 1200);
        io.DeltaTime = 1.0f / 60.0f;
        unsigned char* pixels;
        int width, height;
        io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);

        auto grid = std::make_shared<Grid>();
        grid->setDimensions(1, 1, 1);
        grid->setSpacing(1, 1, 1);
        auto dose = std::make_shared<DoseMatrix>();
        dose->setGrid(*grid);
        dose->allocate();
        dose->at(0, 0, 0) = 60.0;
        state.patientData = std::make_shared<PatientData>();
        auto ct = std::make_unique<Volume<int16_t>>();
        ct->setGrid(*grid);
        ct->allocate();
        state.patientData->setCTVolume(std::move(ct));
        auto structures = std::make_unique<StructureSet>();
        for (int i = 0; i < 12; ++i) {
            auto structure = std::make_unique<Structure>();
            structure->setName("Structure_" + std::to_string(i));
            structure->setType("PTV");
            structure->setVoxelIndices({0});
            structures->addStructure(std::move(structure));
        }
        state.patientData->setStructureSet(std::move(structures));
        state.doseStore.addDose("Dose", dose, grid);
        state.syncSelectedDose();
    }

    void TearDown() override { ImGui::DestroyContext(); }

    void frame(IPanel& panel, ImVec2 size = ImVec2(420, 750)) {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
        ImGui::SetNextWindowSize(size, ImGuiCond_Always);
        panel.render();
        ImGui::Render();
    }

    ImGuiTable* table(const char* windowName, const char* tableName) {
        const auto* window = ImGui::FindWindowByName(windowName);
        if (!window) return nullptr;
        return GImGui->Tables.GetByKey(ImHashStr(tableName, 0, window->ID));
    }
};

TEST_F(PanelLayoutTest, LastDvhStructureWrapsAndCanBeClickedWithoutResizing) {
    DVHPanel panel(state);
    frame(panel);
    frame(panel);
    const auto* window = ImGui::FindWindowByName("DVH");
    ASSERT_NE(window, nullptr);
    const ImGuiID lastToggle = ImHashStr("Structure_11", 0, window->ID);
    auto& io = ImGui::GetIO();
    ImVec2 hit(-1, -1);
    for (int y = 100; y < 500 && hit.x < 0; y += 8) {
        for (int x = 10; x < 400; x += 8) {
            io.AddMousePosEvent(static_cast<float>(x), static_cast<float>(y));
            frame(panel);
            if (GImGui->HoveredId == lastToggle) {
                hit = ImVec2(static_cast<float>(x), static_cast<float>(y));
                break;
            }
        }
    }
    ASSERT_GE(hit.x, 0) << "Last structure toggle must remain inside the narrow panel";
    EXPECT_LE(window->ScrollMax.x, 1.0f);
    io.AddMousePosEvent(-100, -100);
    frame(panel);
    const int withCurve = ImGui::GetDrawData()->TotalVtxCount;
    io.AddMousePosEvent(hit.x, hit.y);
    frame(panel);
    io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
    frame(panel);
    io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
    frame(panel);
    io.AddMousePosEvent(-100, -100);
    frame(panel);
    EXPECT_LT(ImGui::GetDrawData()->TotalVtxCount, withCurve);
}

TEST_F(PanelLayoutTest, StatisticsFitLongNamesCountsAndHeadersAndScrollBothWays) {
    state.patientData->getStructureSet()->getStructure(0)->setName("A_very_long_structure_name_for_statistics");
    state.getOrComputeStats(0, *state.patientData);
    auto& stats = state.statsCache.at(state.doseStore.getSelected()->id);
    stats.front().numVoxels = 1000000;
    stats.front().totalStructureVoxels = 10000000;
    DoseStatsPanel panel(state);
    for (int i = 0; i < 5; ++i) frame(panel, ImVec2(420, 200));
    auto* data = table("Dose Statistics", "DoseStatsTable");
    ASSERT_NE(data, nullptr);
    EXPECT_GE(data->Columns[0].WidthGiven, ImGui::CalcTextSize(stats.front().name.c_str()).x);
    EXPECT_LE(data->Columns[0].WidthGiven, ImGui::CalcTextSize(stats.front().name.c_str()).x + 2.0f);
    EXPECT_GE(data->Columns[1].WidthGiven, ImGui::CalcTextSize("1000000 / 10000000").x);
    EXPECT_GE(data->Columns[2].WidthGiven, ImGui::CalcTextSize("Coverage (%)").x);
    EXPECT_GT(data->InnerWindow->ScrollMax.x, 0);
    EXPECT_GT(data->InnerWindow->ScrollMax.y, 0);
}

TEST_F(PanelLayoutTest, StructureTableFitsCustomAndPlainHeadersWhenTheyAreLongest) {
    ImGui::GetStyle().ItemSpacing.x = 0.0f;
    auto* structures = state.patientData->getStructureSet();
    for (size_t i = 0; i < structures->getCount(); ++i) {
        structures->getStructure(i)->setName("I");
        structures->getStructure(i)->setType("X");
    }

    PatientPanel panel(state, 1, 2);
    for (int i = 0; i < 5; ++i) frame(panel, ImVec2(800, 700));
    auto* data = table("Patient Data", "Structures");
    ASSERT_NE(data, nullptr);
    const float visibilityWidth = ImGui::GetTextLineHeight() * 1.58f +
                                  ImGui::GetStyle().FramePadding.x * 2.0f;
    EXPECT_GE(data->Columns[0].WidthGiven + 1.0f, visibilityWidth + ImGui::CalcTextSize("Name").x);
    EXPECT_GE(data->Columns[1].WidthGiven, ImGui::CalcTextSize("Type").x);
    EXPECT_GE(data->Columns[2].WidthGiven, ImGui::CalcTextSize("Contours").x);

    structures->getStructure(0)->setType("Longer structure type");
    for (int i = 0; i < 5; ++i) frame(panel, ImVec2(800, 700));
    EXPECT_GE(data->Columns[1].WidthGiven, ImGui::CalcTextSize("Longer structure type").x);

    for (int i = 0; i < 5; ++i) frame(panel, ImVec2(140, 700));
    ASSERT_GT(data->InnerWindow->ScrollMax.x, 0.0f);
    ImGui::SetScrollX(data->InnerWindow, data->InnerWindow->ScrollMax.x);
    for (int i = 0; i < 5; ++i) frame(panel, ImVec2(140, 700));
    const auto& contours = data->Columns[2];
    EXPECT_LE(contours.ClipRect.Min.x, contours.WorkMinX);
    EXPECT_GE(contours.ClipRect.Max.x, contours.WorkMinX + ImGui::CalcTextSize("Contours").x);
}

TEST_F(PanelLayoutTest, ObjectiveCombosFitLongestOptionsAndDoseNamesGrowWithData) {
    ImGui::GetIO().FontGlobalScale = 1.5f;
    ImGui::GetStyle().ScaleAllSizes(1.5f);
    const std::string longName = "A_very_long_target_structure_name";
    state.patientData->getStructureSet()->getStructure(0)->setName(longName);
    state.plan = std::make_shared<Plan>();
    state.plan->setDeliveryMode(DeliveryMode::VMAT);
    state.plan->setPrescribedDose(60);
    OptimizationPanel optimization(state);
    for (int i = 0; i < 5; ++i) frame(optimization);
    auto* objectives = table("Optimization", "ObjTable");
    ASSERT_NE(objectives, nullptr);
    const float decoration = ImGui::GetStyle().FramePadding.x * 2 + ImGui::GetFrameHeight();
    EXPECT_GE(objectives->Columns[0].WidthGiven, ImGui::CalcTextSize(longName.c_str()).x + decoration);
    EXPECT_GE(objectives->Columns[1].WidthGiven, ImGui::CalcTextSize("Squared Deviation").x + decoration);
    EXPECT_GT(objectives->InnerWindow->ScrollMax.x, 0);

    PatientPanel patient(state);
    for (int i = 0; i < 5; ++i) frame(patient);
    auto* doses = table("Patient Data", "DoseList");
    ASSERT_NE(doses, nullptr);
    const float initialWidth = doses->Columns[1].WidthGiven;
    const std::string renamed = "Dose_with_a_much_longer_name_than_the_original";
    state.doseStore.renameDose(0, renamed);
    for (int i = 0; i < 5; ++i) frame(patient);
    EXPECT_GT(doses->Columns[1].WidthGiven, initialWidth);
    EXPECT_GE(doses->Columns[1].WidthGiven, ImGui::CalcTextSize(renamed.c_str()).x);
    EXPECT_GT(doses->InnerWindow->ScrollMax.x, 0);
    state.doseStore.renameDose(0, "Dose");
    for (int i = 0; i < 5; ++i) frame(patient);
    EXPECT_FLOAT_EQ(doses->Columns[1].WidthGiven, initialWidth);
}

} // namespace optirad::tests
