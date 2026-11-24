#include <gtest/gtest.h>
#include "gui/panels/PatientPanel.hpp"
#include "gui/panels/DoseStatsPanel.hpp"
#include "gui/panels/DVHPanel.hpp"
#include "io/DicomExporter.hpp"
#include <dcmtk/dcmdata/dctk.h>
#include <imgui.h>
#include <imgui_internal.h>
#include <filesystem>
#include <chrono>
#include <thread>

namespace optirad::tests {

class PatientImportGuiTest : public testing::Test {
protected:
    GuiAppState state;
    std::filesystem::path dir;

    void SetUp() override {
        char uid[100];
        dcmGenerateUniqueIdentifier(uid, SITE_INSTANCE_UID_ROOT);
        dir = std::filesystem::path(testing::TempDir()) / ("optirad_import_" + std::string(uid));
        ASSERT_TRUE(std::filesystem::create_directory(dir));
        writeImage("ct.dcm", UID_CTImageStorage, {0, 0, 0, 0});
        writeImage("dose.dcm", UID_RTDoseStorage, {1000, 2000, 3000, 4000});

        DcmFileFormat file;
        auto* ds = file.getDataset();
        ASSERT_TRUE(ds->putAndInsertString(DCM_SOPClassUID, UID_RTStructureSetStorage).good());
        DcmItem* roi = nullptr;
        ASSERT_TRUE(ds->findOrCreateSequenceItem(DCM_StructureSetROISequence, roi, 0).good());
        ASSERT_TRUE(roi->putAndInsertString(DCM_ROINumber, "1").good());
        ASSERT_TRUE(roi->putAndInsertString(DCM_ROIName, "BODY").good());
        DcmItem* contours = nullptr;
        ASSERT_TRUE(ds->findOrCreateSequenceItem(DCM_ROIContourSequence, contours, 0).good());
        ASSERT_TRUE(contours->putAndInsertString(DCM_ReferencedROINumber, "1").good());
        DcmItem* contour = nullptr;
        ASSERT_TRUE(contours->findOrCreateSequenceItem(DCM_ContourSequence, contour, 0).good());
        ASSERT_TRUE(contour->putAndInsertString(DCM_ContourGeometricType, "CLOSED_PLANAR").good());
        ASSERT_TRUE(contour->putAndInsertString(DCM_NumberOfContourPoints, "4").good());
        ASSERT_TRUE(contour->putAndInsertString(
            DCM_ContourData, "-0.5\\-0.5\\0\\1.5\\-0.5\\0\\1.5\\1.5\\0\\-0.5\\1.5\\0").good());
        ASSERT_TRUE(file.saveFile((dir / "struct.dcm").c_str(), EXS_LittleEndianExplicit).good());

        ImGui::CreateContext();
        auto& io = ImGui::GetIO();
        io.IniFilename = nullptr;
        io.DisplaySize = ImVec2(1600, 1200);
        io.DeltaTime = 1.0f / 60.0f;
        unsigned char* pixels;
        int width, height;
        io.Fonts->GetTexDataAsRGBA32(&pixels, &width, &height);
    }

    void TearDown() override {
        ImGui::DestroyContext();
        std::filesystem::remove_all(dir);
    }

    void writeImage(const char* name, const char* sopClass, const std::vector<Uint16>& pixels) {
        DcmFileFormat file;
        auto* ds = file.getDataset();
        ASSERT_TRUE(ds->putAndInsertString(DCM_SOPClassUID, sopClass).good());
        ASSERT_TRUE(ds->putAndInsertUint16(DCM_Rows, 2).good());
        ASSERT_TRUE(ds->putAndInsertUint16(DCM_Columns, 2).good());
        ASSERT_TRUE(ds->putAndInsertUint16(DCM_BitsAllocated, 16).good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_PixelSpacing, "1\\1").good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_ImagePositionPatient, "0\\0\\0").good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_ImageOrientationPatient, "1\\0\\0\\0\\1\\0").good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_RescaleSlope, "1").good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_RescaleIntercept, "0").good());
        ASSERT_TRUE(ds->putAndInsertString(DCM_DoseGridScaling, "0.01").good());
        ASSERT_TRUE(ds->putAndInsertUint16Array(DCM_PixelData, pixels.data(), pixels.size()).good());
        ASSERT_TRUE(file.saveFile((dir / name).c_str(), EXS_LittleEndianExplicit).good());
    }

    void expectStats(double scale = 1.0) {
        const auto* entry = state.doseStore.getSelected();
        ASSERT_NE(entry, nullptr);
        const auto it = state.statsCache.find(entry->id);
        ASSERT_NE(it, state.statsCache.end());
        ASSERT_EQ(it->second.size(), 1u);
        const auto& stats = it->second.front();
        EXPECT_EQ(stats.name, "BODY");
        EXPECT_EQ(stats.numVoxels, 4u);
        EXPECT_DOUBLE_EQ(stats.minDose, 10.0 * scale);
        EXPECT_DOUBLE_EQ(stats.maxDose, 40.0 * scale);
        EXPECT_DOUBLE_EQ(stats.meanDose, 25.0 * scale);
        EXPECT_DOUBLE_EQ(stats.coveragePercent(), 100.0);
        const auto curves = PlanAnalysis::computeDVHCurves(it->second, 40.0 * scale, 5);
        ASSERT_EQ(curves.size(), 1u);
        EXPECT_FLOAT_EQ(curves.front().volumes.front(), 100.0f);
        EXPECT_FLOAT_EQ(curves.front().volumes.back(), 25.0f);
    }

    bool waitForImport(PatientPanel& panel) {
        const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
        while (panel.isImporting() && std::chrono::steady_clock::now() < deadline) {
            panel.update();
            std::this_thread::sleep_for(std::chrono::milliseconds(1));
        }
        return !panel.isImporting();
    }
};

TEST_F(PatientImportGuiTest, ImportedDoseIsAnalysedAfterBackgroundImportWithoutAPlan) {
    PatientPanel patient(state);
    DoseStatsPanel stats(state);
    DVHPanel dvh(state);
    ImGui::NewFrame();
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(state.dicomLoaded());
    EXPECT_EQ(patient.getPatientData(), state.patientData.get());
    ASSERT_TRUE(state.doseAvailable());
    EXPECT_EQ(state.doseResult, state.doseStore.getSelected()->dose);
    EXPECT_EQ(state.displayGrid, state.doseStore.getSelected()->grid);
    EXPECT_FALSE(state.planCreated());
    EXPECT_FALSE(state.optimizationDone());

    stats.render();
    expectStats();
    state.statsCache.clear();
    dvh.render();
    expectStats();
    ImGui::Render();
}

TEST_F(PatientImportGuiTest, ReimportRefreshesBothPanelsForTheNewDose) {
    PatientPanel patient(state);
    DoseStatsPanel stats(state);
    DVHPanel dvh(state);
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(state.patientData);
    ImGui::NewFrame();
    stats.render();
    dvh.render();
    expectStats();
    ImGui::Render();

    writeImage("dose.dcm", UID_RTDoseStorage, {2000, 4000, 6000, 8000});
    ImGui::NewFrame();
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    EXPECT_EQ(state.doseStore.count(), 1);
    stats.render();
    expectStats(2.0);
    state.statsCache.clear();
    dvh.render();
    expectStats(2.0);
    ImGui::Render();
}

TEST_F(PatientImportGuiTest, ReplacementWithoutDoseClearsOldPatientResults) {
    PatientPanel patient(state);
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(state.patientData);
    auto oldPatient = state.patientData;
    state.plan = std::make_shared<Plan>();
    state.arcCache[1] = {};
    state.getOrComputeStats(0, *state.patientData);
    ASSERT_FALSE(state.statsCache.empty());
    ASSERT_TRUE(std::filesystem::remove(dir / "dose.dcm"));

    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(patient.hasData());
    EXPECT_NE(state.patientData, oldPatient);
    EXPECT_FALSE(state.planCreated());
    EXPECT_FALSE(state.doseAvailable());
    EXPECT_EQ(state.doseStore.count(), 0);
    EXPECT_EQ(state.doseResult, nullptr);
    EXPECT_EQ(state.displayGrid, nullptr);
    EXPECT_TRUE(state.statsCache.empty());
    EXPECT_TRUE(state.arcCache.empty());
}

TEST_F(PatientImportGuiTest, FailedImportPreservesTheCurrentPatientAndDose) {
    PatientPanel patient(state);
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(state.patientData);
    const auto oldPatient = state.patientData;
    const auto oldDose = state.doseResult;
    const int version = state.doseStore.version();

    ASSERT_TRUE(std::filesystem::create_directory(dir / "empty"));
    patient.importDicom((dir / "empty").string());
    ASSERT_TRUE(waitForImport(patient));
    EXPECT_FALSE(state.taskRunning);
    EXPECT_NE(state.taskStatus.find("failed"), std::string::npos);
    EXPECT_EQ(state.patientData, oldPatient);
    EXPECT_EQ(state.doseResult, oldDose);
    EXPECT_EQ(state.doseStore.version(), version);
}

TEST_F(PatientImportGuiTest, ImportCannotReplacePatientWhileATaskIsRunning) {
    PatientPanel patient(state);
    state.taskRunning = true;
    patient.importDicom(dir.string());
    EXPECT_FALSE(patient.hasData());
    EXPECT_FALSE(state.doseAvailable());
}

TEST_F(PatientImportGuiTest, StructureVisibilityButtonsToggleIndividualAndAllStructures) {
    PatientPanel patient(state, 1, 2);
    patient.importDicom(dir.string());
    ASSERT_TRUE(waitForImport(patient));
    ASSERT_TRUE(state.patientData);
    auto* structures = state.patientData->getStructureSet();
    ASSERT_NE(structures, nullptr);
    ASSERT_EQ(structures->getCount(), 1u);
    auto second = std::make_unique<Structure>();
    second->setName("Second");
    structures->addStructure(std::move(second));

    auto frame = [&] {
        ImGui::NewFrame();
        ImGui::SetNextWindowPos(ImVec2(0, 0), ImGuiCond_Always);
        ImGui::SetNextWindowSize(ImVec2(800, 700), ImGuiCond_Always);
        patient.render();
        ImGui::Render();
    };
    frame();
    const auto* window = ImGui::FindWindowByName("Patient Data");
    ASSERT_NE(window, nullptr);
    const ImGuiID tableId = ImHashStr("Structures", 0, window->ID);
    const ImGuiID firstRowId = ImHashStr("##vis0", 0, tableId);
    const ImGuiID allId = ImHashStr("##allVisibility", 0, tableId);

    auto findButton = [&](ImGuiID id) {
        auto& io = ImGui::GetIO();
        for (int y = 30; y < 500; y += 5) {
            for (int x = 10; x < 200; x += 5) {
                io.AddMousePosEvent(static_cast<float>(x), static_cast<float>(y));
                frame();
                if (GImGui->HoveredId == id) return ImVec2(static_cast<float>(x), static_cast<float>(y));
            }
        }
        return ImVec2(-1, -1);
    };
    auto click = [&](ImVec2 pos) {
        auto& io = ImGui::GetIO();
        io.AddMousePosEvent(pos.x, pos.y);
        frame();
        io.AddMouseButtonEvent(ImGuiMouseButton_Left, true);
        frame();
        io.AddMouseButtonEvent(ImGuiMouseButton_Left, false);
        frame();
    };

    const ImVec2 rowButton = findButton(firstRowId);
    ASSERT_GE(rowButton.x, 0);
    click(rowButton);
    EXPECT_FALSE(structures->getStructure(0)->isVisible());
    EXPECT_TRUE(structures->getStructure(1)->isVisible());

    const ImVec2 allButton = findButton(allId);
    ASSERT_GE(allButton.x, 0);
    click(allButton);
    EXPECT_TRUE(structures->getStructure(0)->isVisible());
    EXPECT_TRUE(structures->getStructure(1)->isVisible());

    click(allButton);
    EXPECT_FALSE(structures->getStructure(0)->isVisible());
    EXPECT_FALSE(structures->getStructure(1)->isVisible());

    click(allButton);
    EXPECT_TRUE(structures->getStructure(0)->isVisible());
    EXPECT_TRUE(structures->getStructure(1)->isVisible());
    click(allButton);
    EXPECT_FALSE(structures->getStructure(0)->isVisible());
    EXPECT_FALSE(structures->getStructure(1)->isVisible());

    click(rowButton);
    EXPECT_TRUE(structures->getStructure(0)->isVisible());
    EXPECT_FALSE(structures->getStructure(1)->isVisible());
    click(allButton);
    EXPECT_TRUE(structures->getStructure(0)->isVisible());
    EXPECT_TRUE(structures->getStructure(1)->isVisible());
}

TEST_F(PatientImportGuiTest, PublishesImagesThenStructuresThenPlanAndRefreshesEarlyStats) {
    Plan plan;
    Machine machine;
    MachineMeta meta;
    meta.SAD = 1000.0;
    machine.setMeta(meta);
    plan.setMachine(machine);
    plan.setNumOfFractions(30);
    Stf stf;
    Beam beam;
    beam.setSAD(1000.0);
    beam.setIsocenter({0.0, 0.0, 0.0});
    stf.addBeam(beam);
    LeafSequenceResult sequence;
    sequence.leafPairBoundariesZ = {-5.0, 5.0};
    Aperture segment;
    segment.bankA = {-10.0};
    segment.bankB = {10.0};
    segment.weight = 10.0;
    sequence.segments.push_back(segment);
    sequence.totalMU = 10.0;
    ASSERT_TRUE(DicomExporter().exportRTPlan(plan, stf, {sequence}, DicomContext{}, dir.string()));

    PatientPanel patient(state);
    ImGui::NewFrame();
    ImGui::SetNextWindowSize(ImVec2(800, 700), ImGuiCond_Always);
    patient.render();
    ImGui::Render();
    patient.importDicom(dir.string());
    ASSERT_TRUE(patient.isImporting());
    ASSERT_TRUE(state.taskRunning);
    EXPECT_FALSE(state.patientData); // worker results are only published by update()
    EXPECT_FALSE(state.taskStatus.empty());
    patient.importDicom((dir / "empty").string()); // a second import cannot change the active one

    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    bool sawImages = false;
    bool sawStructures = false;
    while (patient.isImporting() && std::chrono::steady_clock::now() < deadline) {
        ImGui::NewFrame();
        patient.render();
        ImGui::Render();
        EXPECT_GT(ImGui::GetDrawData()->TotalVtxCount, 0);
        if (state.patientData && !state.patientData->hasStructures()) {
            sawImages = true;
            ASSERT_TRUE(state.doseAvailable());
            EXPECT_TRUE(state.patientData->hasValidCT());
            EXPECT_FALSE(state.patientData->hasImportedPlan());
            EXPECT_TRUE(state.taskRunning);
            EXPECT_TRUE(state.getOrComputeStats(0, *state.patientData).empty());
        }
        if (state.patientData && state.patientData->hasStructures() && !state.patientData->hasImportedPlan()) {
            sawStructures = true;
            EXPECT_TRUE(sawImages);
            EXPECT_TRUE(state.taskRunning);
            state.getOrComputeStats(0, *state.patientData);
            expectStats();
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    ASSERT_FALSE(patient.isImporting());
    EXPECT_TRUE(sawImages);
    EXPECT_TRUE(sawStructures);
    EXPECT_FALSE(state.taskRunning);
    ASSERT_TRUE(state.patientData->hasImportedPlan());
    ASSERT_TRUE(state.importedPlan);
    EXPECT_EQ(state.plan, state.importedPlan);
    EXPECT_EQ(state.dosePlans.at(state.doseStore.getSelected()->id), state.importedPlan);
    EXPECT_FLOAT_EQ(state.taskProgress, 1.0f);
}

TEST_F(PatientImportGuiTest, HiddenPanelStillFinishesImportAndReleasesTaskGuard) {
    PatientPanel patient(state);
    patient.setVisible(false);
    patient.importDicom(dir.string());
    const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(10);
    while (patient.isImporting() && std::chrono::steady_clock::now() < deadline) {
        ImGui::NewFrame();
        patient.render();
        ImGui::Render();
        std::this_thread::sleep_for(std::chrono::milliseconds(1));
    }
    EXPECT_FALSE(patient.isImporting());
    EXPECT_FALSE(state.taskRunning);
    ASSERT_TRUE(state.patientData);
    EXPECT_TRUE(state.patientData->hasStructures());
}

TEST_F(PatientImportGuiTest, DestroyingAnImportingPanelJoinsItsWorkerAndReleasesTaskGuard) {
    {
        PatientPanel patient(state);
        patient.importDicom(dir.string());
        ASSERT_TRUE(patient.isImporting());
        ASSERT_TRUE(state.taskRunning);
    }
    EXPECT_FALSE(state.taskRunning);
}

TEST_F(PatientImportGuiTest, SynchronousImporterRetainsCompleteImportAndClearsDoseOnReuse) {
    DicomImporter importer;
    auto data = importer.importAll(dir.string());
    ASSERT_TRUE(data);
    EXPECT_TRUE(data->hasValidCT());
    EXPECT_TRUE(data->hasStructures());
    EXPECT_TRUE(data->getImportedDoseGridResolution().has_value());
    ASSERT_TRUE(std::filesystem::remove(dir / "dose.dcm"));
    data = importer.importAll(dir.string());
    ASSERT_TRUE(data);
    EXPECT_TRUE(data->hasStructures());
    EXPECT_FALSE(data->getImportedDoseGridResolution().has_value());
    const auto [dose, grid] = importer.takeImportedDose();
    EXPECT_FALSE(dose);
    EXPECT_FALSE(grid);
}

} // namespace optirad::tests
