#include <gtest/gtest.h>
#include "gui/views/BevView.hpp"
#include "gui/views/renderers/BevRenderer.hpp"
#include "gui/views/SliceView.hpp"
#include "gui/views/View3D.hpp"
#include "gui/views/renderers/BeamRenderer.hpp"
#include "gui/panels/DeliveryPanel.hpp"
#include "sequencing/ArcSequencer.hpp"
#include <GLFW/glfw3.h>
#include <imgui_internal.h>
#include <array>
#include <cmath>
#include <memory>

namespace optirad::tests {
namespace {

class BevRendererTest : public testing::Test {
protected:
    GLFWwindow* window = nullptr;
    Beam beam;
    Aperture aperture;
    const std::vector<double> boundaries{-30.0, 30.0};

    void SetUp() override {
        if (!glfwInit()) GTEST_SKIP() << "OpenGL display unavailable";
        glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
        glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
        window = glfwCreateWindow(512, 512, "BEV test", nullptr, nullptr);
        if (!window) {
            glfwTerminate();
            GTEST_SKIP() << "OpenGL 4.0 context unavailable";
        }
        glfwMakeContextCurrent(window);
        glewExperimental = GL_TRUE;
        ASSERT_EQ(glewInit(), GLEW_OK);
        while (glGetError() != GL_NO_ERROR) {}
        beam.setSAD(1000.0);
        beam.setSCD(500.0);
        aperture.bankA = {-20.0};
        aperture.bankB = {20.0};
    }

    void TearDown() override {
        if (ImGui::GetCurrentContext()) ImGui::DestroyContext();
        if (window) {
            glfwDestroyWindow(window);
            glfwTerminate();
        }
    }

    GLuint render(BevRenderer& renderer, PatientData* patient = nullptr) {
        return renderer.render(patient, beam, {&aperture}, boundaries, 100.0, 100.0,
                               512, 512, 1.0f);
    }

    std::array<unsigned char, 4> pixel(GLuint texture, int x, int y) {
        std::vector<unsigned char> image(512 * 512 * 4);
        glBindTexture(GL_TEXTURE_2D, texture);
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
        const auto offset = static_cast<size_t>(y * 512 + x) * 4;
        return {image[offset], image[offset + 1], image[offset + 2], image[offset + 3]};
    }

    static std::shared_ptr<PatientData> patient() {
        auto data = std::make_shared<PatientData>();
        auto structures = std::make_unique<StructureSet>();
        auto target = std::make_unique<Structure>();
        target->setName("PTV");
        target->setColor(255, 40, 20);
        for (double z = -20.0; z <= 20.0; z += 1.0) {
            Contour contour;
            contour.zPosition = z;
            contour.points = {{-15.0, -15.0, z}, {15.0, -15.0, z},
                              {15.0, 15.0, z}, {-15.0, 15.0, z}};
            target->addContour(contour);
        }
        structures->addStructure(std::move(target));
        data->setStructureSet(std::move(structures));
        return data;
    }

    void deliveryState(GuiAppState& state) {
        state.patientData = patient();
        Grid grid;
        grid.setDimensions(64, 64, 64);
        grid.setSpacing(4.0, 4.0, 4.0);
        grid.setOrigin({-128.0, -128.0, -128.0});
        auto ct = std::make_unique<Volume<int16_t>>();
        ct->setGrid(grid);
        ct->allocate();
        state.patientData->setCTVolume(std::move(ct));
        state.stf = std::make_shared<Stf>();
        beam.generateRays(5.0, {60.0, 60.0});
        state.stf->addBeam(beam);
        auto arc = std::make_shared<ArcPipelineResult>();
        arc->controlPointStf = state.stf;
        DeliveredControlPoint cp;
        cp.geometry.gantryAngle = 180.0;
        cp.geometry.width = 360.0;
        cp.timeSec = 1.0;
        cp.mu = 1.0;
        cp.gantrySpeed = 360.0;
        cp.aperture = aperture;
        cp.leafPairBoundariesZ = boundaries;
        arc->sequence.controlPoints = {cp};
        arc->sequence.totalTimeSec = 1.0;
        const int id = state.doseStore.addDose("Delivery", nullptr, nullptr);
        state.arcCache[id].result = arc;
        state.delivery.synchronize(state);
    }

    static bool yellow(ImU32 color) {
        return ((color >> IM_COL32_R_SHIFT) & 255) > 240 &&
               ((color >> IM_COL32_G_SHIFT) & 255) > 200 &&
               ((color >> IM_COL32_B_SHIFT) & 255) < 100;
    }
};

TEST_F(BevRendererTest, SourceProjectionIsIndependentOfCollimatorDistanceAndBeamRotation) {
    BevRenderer renderer;
    for (double scd : {300.0, 500.0, 700.0}) {
        beam.setSCD(scd);
        beam.setIsocenter(121.0, -55.0, 260.0);
        beam.setGantryAngle(scd / 2.0);
        beam.setCouchAngle(35.0);
        const auto texture = render(renderer);
        ASSERT_NE(texture, 0u);
        EXPECT_LT(pixel(texture, 256, 256)[0], 20); // aperture
        EXPECT_GT(pixel(texture, 380, 256)[0], 50); // leaf
        const int tipX = static_cast<int>(256.0 + 20.0 / 115.0 * 256.0);
        EXPECT_GT(pixel(texture, tipX, 256)[1], 200); // projected right tip
        EXPECT_EQ(glGetError(), GL_NO_ERROR);
    }
}

TEST_F(BevRendererTest, ZoomKeepsSourcePerspectiveAlignmentAtEveryMlcDistance) {
    BevRenderer renderer;
    for (double scd : {300.0, 500.0, 700.0}) {
        beam.setSCD(scd);
        const auto texture = renderer.render(nullptr, beam, {&aperture}, boundaries,
            100.0, 100.0, 512, 512, 2.0f);
        const int tipX = static_cast<int>(256.0 + 20.0 / (115.0 / 2.0) * 256.0);
        EXPECT_GT(pixel(texture, tipX, 256)[1], 200);
        EXPECT_LT(pixel(texture, 256, 256)[0], 20);
        EXPECT_EQ(glGetError(), GL_NO_ERROR);
    }
}

TEST_F(BevRendererTest, FirstFrameAndResizeRenderAnatomyIntoTheRequestedFramebuffer) {
    BevRenderer renderer;
    auto data = patient();
    glViewport(11, 12, 123, 234);
    glEnable(GL_SCISSOR_TEST);
    glScissor(0, 0, 1, 1);
    glEnable(GL_BLEND);
    glDisable(GL_DEPTH_TEST);
    const auto texture = render(renderer, data.get());
    const auto visible = pixel(texture, 256, 256);
    EXPECT_GT(visible[0], visible[1] + 30);
    GLint fbo, viewport[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &fbo);
    glGetIntegerv(GL_VIEWPORT, viewport);
    EXPECT_EQ(fbo, 0);
    EXPECT_EQ(viewport[0], 11);
    EXPECT_EQ(viewport[2], 123);
    EXPECT_TRUE(glIsEnabled(GL_SCISSOR_TEST));
    EXPECT_TRUE(glIsEnabled(GL_BLEND));
    EXPECT_FALSE(glIsEnabled(GL_DEPTH_TEST));
    data->getStructureSet()->getStructure(0)->setVisible(false);
    EXPECT_LT(pixel(render(renderer, data.get()), 256, 256)[0], 20);
    data->getStructureSet()->getStructure(0)->setVisible(true);
    renderer.render(data.get(), beam, {&aperture}, boundaries, 100.0, 100.0,
                    256, 300, 1.0f);
    EXPECT_GT(pixel(render(renderer, data.get()), 256, 256)[0], 50);
    EXPECT_EQ(glGetError(), GL_NO_ERROR);
}

TEST_F(BevRendererTest, CombinedSegmentsUseOutlinesInsteadOfOverlappingLeafBanks) {
    BevRenderer renderer;
    Aperture second = aperture;
    second.bankA[0] = -10.0;
    second.bankB[0] = 10.0;
    const auto texture = renderer.render(nullptr, beam, {&aperture, &second}, boundaries,
        100.0, 100.0, 512, 512, 1.0f);
    EXPECT_LT(pixel(texture, 380, 256)[0], 20);
    EXPECT_GT(pixel(texture, static_cast<int>(256.0 + 10.0 / 115.0 * 256.0), 256)[1], 200);
    EXPECT_EQ(glGetError(), GL_NO_ERROR);
}

TEST_F(BevRendererTest, SmlcAndCachedVmatSelectionsRenderWithoutOptimizationWeights) {
    ImGui::CreateContext();
    auto& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.DisplaySize = ImVec2(1000, 1000);
    io.DeltaTime = 1.0f / 60.0f;
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    GuiAppState state;
    state.patientData = patient();
    state.stf = std::make_shared<Stf>();
    beam.generateRays(5.0, {60.0, 60.0});
    state.stf->addBeam(beam);
    LeafSequenceResult sequence;
    sequence.segments = {aperture};
    sequence.leafPairBoundariesZ = boundaries;
    state.leafSequences = {sequence};
    const int smlc = state.doseStore.addDose("SMLC", nullptr, nullptr);
    state.seqCache[smlc].sequences = state.leafSequences;
    {
        BevView view(state);
        auto frame = [&] {
            ImGui::NewFrame();
            ImGui::LogToBuffer();
            ImGui::SetNextWindowSize(ImVec2(900, 900));
            view.render();
            const std::string text = GImGui->LogBuffer.c_str();
            ImGui::LogFinish();
            ImGui::Render();
            EXPECT_GT(ImGui::GetDrawData()->TotalVtxCount, 0);
            EXPECT_EQ(glGetError(), GL_NO_ERROR);
            return text;
        };
        frame();
        frame();
        auto arc = std::make_shared<ArcPipelineResult>();
        arc->controlPointStf = state.stf;
        arc->sequence.controlPoints.resize(1);
        arc->sequence.controlPoints[0].aperture = aperture;
        arc->sequence.controlPoints[0].leafPairBoundariesZ = boundaries;
        arc->sequence.controlPoints[0].geometry.width = 2.0;
        arc->sequence.controlPoints[0].mu = 1.0;
        ArcSequencer::computeTiming(arc->sequence, MachineConstraints{});
        const int vmat = state.doseStore.addDose("VMAT", nullptr, nullptr);
        state.arcCache[vmat].result = arc;
        state.arcResult.reset(); // the selected cache, not the last optimization, is authoritative
        const auto fullSpeed = frame();
        EXPECT_NE(fullSpeed.find("Gantry speed: 6.00 deg/s (delivery)"), std::string::npos);
        EXPECT_NE(fullSpeed.find("Playback"), std::string::npos);
        arc->sequence.controlPoints[0].mu = 40.0;
        ArcSequencer::computeTiming(arc->sequence, MachineConstraints{});
        EXPECT_NE(frame().find("Gantry speed: 0.50 deg/s (delivery)"), std::string::npos);
        state.doseStore.setSelected(0);
        frame();
    }
    ImGui::DestroyContext();
}

TEST_F(BevRendererTest, SceneMouseDragsDoNotChangeProjectionButScrollingZooms) {
    ImGui::CreateContext();
    auto& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.DisplaySize = ImVec2(1000, 1000);
    io.DeltaTime = 1.0f / 60.0f;
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    GuiAppState state;
    state.patientData = patient();
    state.stf = std::make_shared<Stf>();
    beam.generateRays(5.0, {60.0, 60.0});
    state.stf->addBeam(beam);
    {
        BevView view(state);
        auto frame = [&]() {
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(ImVec2(900, 900));
            view.render();
            ImGui::Render();
            EXPECT_EQ(glGetError(), GL_NO_ERROR);
            std::vector<unsigned char> image;
            const auto* drawData = ImGui::GetDrawData();
            for (int list = 0; list < drawData->CmdListsCount; ++list) {
                for (const auto& command : drawData->CmdLists[list]->CmdBuffer) {
                    const auto texture = static_cast<GLuint>(command.GetTexID());
                    if (!glIsTexture(texture)) continue;
                    glBindTexture(GL_TEXTURE_2D, texture);
                    GLint w, h;
                    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_WIDTH, &w);
                    glGetTexLevelParameteriv(GL_TEXTURE_2D, 0, GL_TEXTURE_HEIGHT, &h);
                    image.resize(static_cast<size_t>(w) * h * 4);
                    glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, image.data());
                }
            }
            EXPECT_FALSE(image.empty());
            EXPECT_EQ(glGetError(), GL_NO_ERROR);
            return image;
        };
        frame();
        io.AddMousePosEvent(450.0f, 450.0f);
        frame();
        const auto baseline = frame();
        for (int button : {ImGuiMouseButton_Left, ImGuiMouseButton_Middle, ImGuiMouseButton_Right}) {
            io.AddMouseButtonEvent(button, true);
            frame();
            io.AddMousePosEvent(550.0f, 500.0f);
            EXPECT_EQ(frame(), baseline);
            io.AddMouseButtonEvent(button, false);
            frame();
            io.AddMousePosEvent(450.0f, 450.0f);
            frame();
        }
        io.AddMouseWheelEvent(0.0f, 1.0f);
        EXPECT_NE(frame(), baseline);
    }
    ImGui::DestroyContext();
}

TEST_F(BevRendererTest, DeliveryProjectionIsVisibleInEverySliceAndHiddenByOneCheckboxState) {
    ImGui::CreateContext();
    auto& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.DisplaySize = ImVec2(1000, 1000);
    io.DeltaTime = 1.0f / 60.0f;
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    GuiAppState state;
    deliveryState(state);
    state.delivery.seek(0.125);
    for (auto orientation : {SliceOrientation::Axial, SliceOrientation::Coronal, SliceOrientation::Sagittal}) {
        SliceView view(orientation);
        view.setPatientData(state.patientData.get());
        view.setDelivery(&state.delivery);
        auto frame = [&] {
            ImGui::NewFrame();
            ImGui::SetNextWindowPos(ImVec2(0, 0));
            ImGui::SetNextWindowSize(ImVec2(900, 900));
            view.render();
            ImGui::Render();
            std::vector<std::array<float, 2>> yellowVertices;
            const auto* data = ImGui::GetDrawData();
            for (int i = 0; i < data->CmdListsCount; ++i) {
                for (const auto& vertex : data->CmdLists[i]->VtxBuffer) {
                    EXPECT_TRUE(std::isfinite(vertex.pos.x));
                    EXPECT_TRUE(std::isfinite(vertex.pos.y));
                    if (yellow(vertex.col)) yellowVertices.push_back({vertex.pos.x, vertex.pos.y});
                }
            }
            EXPECT_EQ(glGetError(), GL_NO_ERROR);
            return yellowVertices;
        };
        state.delivery.showOverlays = true;
        frame();
        const auto initial = frame();
        EXPECT_GT(initial.size(), 400u);
        state.delivery.seek(0.25);
        EXPECT_NE(frame(), initial);
        state.delivery.seek(0.125);
        state.delivery.showOverlays = false;
        EXPECT_TRUE(frame().empty());
        EXPECT_DOUBLE_EQ(state.delivery.time(), 0.125);
    }
}

TEST_F(BevRendererTest, DeliveryPanelAndBevReadOneClockWithoutAdvancingIt) {
    ImGui::CreateContext();
    auto& io = ImGui::GetIO();
    io.IniFilename = nullptr;
    io.DisplaySize = ImVec2(1600, 1000);
    io.DeltaTime = 1.0f / 60.0f;
    unsigned char* pixels;
    int width, height;
    io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    GuiAppState state;
    deliveryState(state);
    DeliveryPanel panel(state);
    BevView bev(state);
    state.delivery.togglePlaying();
    state.delivery.advance(0.25);
    for (int i = 0; i < 2; ++i) {
        ImGui::NewFrame();
        ImGui::LogToBuffer();
        ImGui::SetNextWindowPos(ImVec2(0, 0));
        ImGui::SetNextWindowSize(ImVec2(600, 900));
        panel.render();
        ImGui::SetNextWindowPos(ImVec2(650, 0));
        ImGui::SetNextWindowSize(ImVec2(900, 900));
        bev.render();
        const std::string text = GImGui->LogBuffer.c_str();
        ImGui::LogFinish();
        ImGui::Render();
        const auto first = text.find("Gantry 90.0");
        ASSERT_NE(first, std::string::npos);
        EXPECT_NE(text.find("Gantry 90.0", first + 1), std::string::npos);
        EXPECT_DOUBLE_EQ(state.delivery.time(), 0.25);
        EXPECT_TRUE(state.delivery.playing);
        EXPECT_EQ(glGetError(), GL_NO_ERROR);
    }
}

TEST_F(BevRendererTest, View3DShowsMovingYellowDeliveryAndStfBeamsStartHidden) {
    GuiAppState state;
    deliveryState(state);
    View3D view;
    view.init();
    view.setPatientData(state.patientData.get());
    view.setStf(state.stf.get());
    view.setDelivery(&state.delivery);
    auto image = [&] {
        view.render(512, 512);
        std::vector<unsigned char> pixels(512 * 512 * 4);
        glBindTexture(GL_TEXTURE_2D, view.getTextureID());
        glGetTexImage(GL_TEXTURE_2D, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
        EXPECT_EQ(glGetError(), GL_NO_ERROR);
        return pixels;
    };
    auto countYellow = [](const std::vector<unsigned char>& pixels) {
        size_t count = 0;
        for (size_t i = 0; i < pixels.size(); i += 4)
            if (pixels[i] > 120 && pixels[i + 1] > 100 && pixels[i + 2] < 100) ++count;
        return count;
    };
    const auto visible = image();
    EXPECT_EQ(view.getBeamRenderer()->getBeamCount(), 1u);
    EXPECT_FALSE(view.getBeamRenderer()->isBeamVisible(0));
    EXPECT_GT(countYellow(visible), 100u);
    state.delivery.seek(0.25);
    EXPECT_NE(image(), visible);
    state.delivery.showOverlays = false;
    const auto hidden = image();
    EXPECT_LT(countYellow(hidden), countYellow(visible));
    view.getBeamRenderer()->setBeamVisible(0, true);
    EXPECT_TRUE(view.getBeamRenderer()->isBeamVisible(0));
    EXPECT_NE(image(), hidden);
    view.cleanup();
}

} // namespace
} // namespace optirad::tests
