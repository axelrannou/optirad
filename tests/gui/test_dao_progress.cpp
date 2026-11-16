#include <gtest/gtest.h>
#include "gui/panels/DaoOptimizationView.hpp"
#include <imgui.h>
#include <cmath>
#include <limits>

namespace optirad::tests {

class DaoProgressGuiTest : public testing::Test {
protected:
    GuiAppState state;

    void SetUp() override {
        ImGui::CreateContext();
        auto& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(1200, 2200);
        io.DeltaTime = 1.0f / 60.0f;
        unsigned char* pixels;
        int width, height;
        io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
        state.plan = std::make_shared<Plan>();
        state.plan->setDeliveryMode(DeliveryMode::VMAT);
        state.arcResult = std::make_shared<ArcPipelineResult>();
    }

    void TearDown() override { ImGui::DestroyContext(); }

    int render(DaoOptimizationView& panel) {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(1100, 2100));
        ImGui::Begin("DAO");
        panel.renderProgress();
        ImGui::End();
        ImGui::Render();
        const auto* data = ImGui::GetDrawData();
        for (int i = 0; i < data->CmdListsCount; ++i) {
            for (const auto& vertex : data->CmdLists[i]->VtxBuffer) {
                EXPECT_TRUE(std::isfinite(vertex.pos.x));
                EXPECT_TRUE(std::isfinite(vertex.pos.y));
            }
        }
        return data->TotalVtxCount;
    }
};

TEST_F(DaoProgressGuiTest, DrawsChartsAndBarsForTheSelectedArcHistory) {
    DaoOptimizationView panel(state);
    render(panel);
    const int withoutHistory = render(panel);
    ArcProgressSnapshot sample;
    sample.maxIterations = 50;
    sample.controlPoints = 9;
    sample.progress.totalMu = 200.0;
    sample.breakdown.total = 100.0;
    sample.breakdown.maxDose = 60.0;
    sample.breakdown.contributions = {{"PTV / MinDVH", 80.0}, {"NTO (start 105%)", 20.0}};
    state.arcResult->progressHistory.push_back(sample);
    sample.progress.iteration = 10;
    sample.progress.acceptedSteps = 10;
    sample.progress.totalMu = 190.0;
    sample.breakdown.total = 50.0;
    sample.breakdown.contributions = {{"PTV / MinDVH", 40.0}, {"NTO (start 105%)", 10.0}};
    state.arcResult->progressHistory.push_back(sample);
    EXPECT_GT(render(panel), withoutHistory);
    EXPECT_DOUBLE_EQ(state.arcResult->progressHistory.back().breakdown.total, 50.0);
    state.arcResult.reset();
    EXPECT_LE(render(panel), withoutHistory);
}

TEST_F(DaoProgressGuiTest, ZeroAndInvalidCostsDoNotProduceInvalidGeometry) {
    DaoOptimizationView panel(state);
    ArcProgressSnapshot sample;
    sample.maxIterations = 1;
    sample.breakdown.contributions = {{"PTV", 0.0}};
    state.arcResult->progressHistory = {sample};
    render(panel);
    EXPECT_GT(render(panel), 0);
    sample.breakdown.total = std::numeric_limits<double>::quiet_NaN();
    state.arcResult->progressHistory = {sample};
    EXPECT_GT(render(panel), 0);
    sample.breakdown.total = 1.0;
    sample.breakdown.contributions.front().second = std::numeric_limits<double>::quiet_NaN();
    state.arcResult->progressHistory = {sample};
    EXPECT_GT(render(panel), 0);
}

} // namespace optirad::tests
