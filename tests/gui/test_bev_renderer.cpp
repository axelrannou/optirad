#include <gtest/gtest.h>
#include "gui/views/BevView.hpp"
#include "gui/views/renderers/BevRenderer.hpp"
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

} // namespace
} // namespace optirad::tests
