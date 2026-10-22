#include <gtest/gtest.h>

#include "core/workflow/PlanBuilder.hpp"
#include "dose/FluenceDoseOperator.hpp"
#include "dose/engines/PencilBeamEngine.hpp"
#include "validation/WaterPhantom.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/ObjectiveBuilder.hpp"
#include "dose/CudaRayTracer.hpp"
#include "dose/SiddonRayTracer.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>

namespace optirad::tests {

namespace {

// Three beams on a water cube, shared by the tests.
struct ModelFixture {
    std::shared_ptr<PatientData> patient = WaterPhantom::buildCentredPatient(5.0, 80.0, 20.0);
    PlanBuilder::Result built;
    std::unique_ptr<Grid> grid;
    PencilBeamEngine engine;
    std::vector<double> weights;
    std::vector<size_t> offsets;

    ModelFixture() {
        PlanConfig pc;
        pc.bixelWidth = 5.0;
        pc.gantryStart = 0.0;
        pc.gantryStep = 120.0;
        pc.gantryStop = 360.0;
        built = PlanBuilder::build(pc, patient);
        grid = std::make_unique<Grid>(Grid::createDoseGrid(patient->getGrid(), {5.0, 5.0, 5.0}));

        DoseCalcOptions o;
        o.relativeThreshold = 0.0;
        o.absoluteThreshold = 0.0;
        o.excludeExternal = false;
        o.divergentBixels = true;
        engine.setOptions(o);

        offsets.push_back(0);
        for (size_t b = 0; b < built.stf->getCount(); ++b)
            offsets.push_back(offsets.back() + built.stf->getBeam(b)->getTotalNumOfBixels());
        weights.resize(offsets.back());
        for (size_t i = 0; i < weights.size(); ++i) weights[i] = 1.0 + 0.5 * std::sin(0.37 * static_cast<double>(i));
    }

    std::vector<std::vector<float>> images(const FluenceDoseModel& model) const {
        std::vector<std::vector<float>> out(built.stf->getCount());
        for (size_t b = 0; b < out.size(); ++b)
            model.rasterizeBixels(*built.stf->getBeam(b), weights.data() + offsets[b], out[b]);
        return out;
    }
};

} // namespace

// The matrix-free dose must reproduce the stored-Dij dose of the same bixel weights; the residual is the
// pixel discretisation of the fluence image.
TEST(FluenceDoseModelTest, MatchesTheDijDose) {
    ModelFixture f;
    const auto dij = f.engine.calculateDij(*f.built.plan, *f.built.stf, *f.patient, *f.grid);
    const auto ref = dij.computeDose(f.weights);
    const double maxRef = *std::max_element(ref.begin(), ref.end());
    ASSERT_GT(maxRef, 0.0);

    struct Case { double pixelMm, maxErr, meanErr; };
    for (const Case c : {Case{0.5, 0.01, 0.003}, Case{1.5, 0.04, 0.015}}) {
        FluenceModelOptions fo;
        fo.pixelMm = c.pixelMm;
        const auto model = f.engine.buildFluenceModel(*f.built.plan, *f.built.stf, *f.patient, *f.grid, fo);
        const auto dose = model.computeDose(f.images(model));
        ASSERT_EQ(dose.size(), ref.size());

        double maxErr = 0.0, sumRel = 0.0;
        size_t count = 0;
        for (size_t v = 0; v < ref.size(); ++v) {
            maxErr = std::max(maxErr, std::abs(dose[v] - ref[v]));
            if (ref[v] > 0.2 * maxRef) { sumRel += std::abs(dose[v] - ref[v]) / ref[v]; ++count; }
        }
        ASSERT_GT(count, 100u);
        const double meanRel = sumRel / static_cast<double>(count);
        std::cout << "  pixel " << c.pixelMm << " mm: max error " << 100.0 * maxErr / maxRef
                  << " % of max, mean relative error above 20 % of max " << 100.0 * meanRel << " %, model "
                  << model.memoryBytes() / 1024 << " kB\n";
        EXPECT_LT(maxErr / maxRef, c.maxErr) << "pixel " << c.pixelMm;
        EXPECT_LT(meanRel, c.meanErr) << "pixel " << c.pixelMm;
    }
}

// <A f, g> == <f, A^T g> for the forward dose and its gradient.
TEST(FluenceDoseModelTest, GradientIsTheAdjointOfTheDose) {
    ModelFixture f;
    const auto model = f.engine.buildFluenceModel(*f.built.plan, *f.built.stf, *f.patient, *f.grid);
    const auto fluence = f.images(model);
    const auto dose = model.computeDose(fluence);

    std::vector<double> g(model.numVoxels());
    for (size_t v = 0; v < g.size(); ++v) g[v] = std::cos(0.01 * static_cast<double>(v));
    const auto grads = model.computeGradient(g);

    double lhs = 0.0, rhs = 0.0;
    for (size_t v = 0; v < g.size(); ++v) lhs += dose[v] * g[v];
    for (size_t b = 0; b < fluence.size(); ++b)
        for (size_t p = 0; p < fluence[b].size(); ++p) rhs += static_cast<double>(fluence[b][p]) * grads[b][p];
    EXPECT_NEAR(lhs, rhs, 1e-4 * std::max(1.0, std::abs(lhs)));
}

TEST(FluenceDoseModelTest, ZeroFluenceGivesZeroDoseAndTimingIsReported) {
    ModelFixture f;
    const auto model = f.engine.buildFluenceModel(*f.built.plan, *f.built.stf, *f.patient, *f.grid);
    std::vector<std::vector<float>> empty(model.numBeams());
    const auto zero = model.computeDose(empty);
    EXPECT_EQ(*std::max_element(zero.begin(), zero.end()), 0.0);

    const auto fluence = f.images(model);
    const auto t0 = std::chrono::steady_clock::now();
    for (int i = 0; i < 10; ++i) (void)model.computeDose(fluence);
    const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count() / 10.0;
    std::cout << "  forward dose " << ms << " ms for " << model.numBeams() << " beams\n";
}

// The operator maps bixel weights to dose and back: <A w, g> == <w, A^T g>.
TEST(FluenceDoseOperatorTest, TransposeIsTheAdjointOfTheBixelDose) {
    ModelFixture f;
    auto model = std::make_shared<FluenceDoseModel>(
        f.engine.buildFluenceModel(*f.built.plan, *f.built.stf, *f.patient, *f.grid));
    FluenceDoseOperator op(model, std::make_shared<Stf>(*f.built.stf));
    ASSERT_EQ(op.getNumBixels(), f.weights.size());

    const auto dose = op.computeDose(f.weights);
    std::vector<double> g(op.getNumVoxels());
    for (size_t v = 0; v < g.size(); ++v) g[v] = std::cos(0.013 * static_cast<double>(v));
    std::vector<double> grad(op.getNumBixels(), 0.0);
    op.accumulateTransposeProduct(g, grad);

    double lhs = 0.0, rhs = 0.0;
    for (size_t v = 0; v < g.size(); ++v) lhs += dose[v] * g[v];
    for (size_t i = 0; i < grad.size(); ++i) rhs += f.weights[i] * grad[i];
    EXPECT_NEAR(lhs, rhs, 1e-4 * std::max(1.0, std::abs(lhs)));
    EXPECT_GT(op.getMaxValue(), 0.0);
    EXPECT_TRUE(op.engineOptions().divergentBixels);
}

TEST(CudaFluenceTest, ForwardAndAdjointMatchCpuWithAsymmetricKernels) {
#ifndef OPTIRAD_HAS_CUDA
    GTEST_SKIP() << "CUDA backend not compiled";
#else
    FluenceDoseModel cpu(6, 8, 1.0, 1.0);
    std::array<std::vector<double>, 3> kernel;
    for (size_t component = 0; component < 3; ++component) {
        kernel[component].resize(9);
        for (size_t index = 0; index < 9; ++index)
            kernel[component][index] = (index + 1.0) * (component + 1.0) / 31.0;
    }
    cpu.addKernel(0, kernel, 3);
    for (auto& table : kernel) std::reverse(table.begin(), table.end());
    cpu.addKernel(1, kernel, 3);
    std::vector<FluenceDoseModel::VoxelRecord> records;
    for (uint32_t voxel = 0; voxel < 6; ++voxel)
        records.push_back({voxel, 0.25f + voxel, 6.75f - voxel, {0.7f, 0.2f, 0.1f}});
    cpu.addBeam(0, records);
    cpu.addBeam(1, records);
    cpu.addBeam(0, {});
    auto gpu = cpu;
    ASSERT_TRUE(gpu.enableCuda(0, true));
    EXPECT_TRUE(gpu.cudaGeometryResident());
    std::vector<std::vector<float>> images(3, std::vector<float>(64));
    for (size_t beam = 0; beam < images.size(); ++beam)
        for (size_t pixel = 0; pixel < 64; ++pixel)
            images[beam][pixel] = static_cast<float>(std::sin(0.17 * pixel + beam));
    const std::vector<double> gradient{1.0, -0.2, 0.7, -0.5, 0.1, 0.4};
    const auto expectedDose = cpu.computeDose(images);
    const auto expectedGradient = cpu.computeGradient(gradient);
    auto verify = [&](const FluenceDoseModel& model) {
        const auto dose = model.computeDose(images);
        const auto adjoint = model.computeGradient(gradient);
        double forwardDot = 0.0, adjointDot = 0.0;
        for (size_t voxel = 0; voxel < dose.size(); ++voxel) {
            EXPECT_NEAR(dose[voxel], expectedDose[voxel], 2e-6);
            forwardDot += dose[voxel] * gradient[voxel];
        }
        for (size_t beam = 0; beam < images.size(); ++beam)
            for (size_t pixel = 0; pixel < 64; ++pixel) {
                EXPECT_NEAR(adjoint[beam][pixel], expectedGradient[beam][pixel], 2e-6);
                adjointDot += images[beam][pixel] * adjoint[beam][pixel];
            }
        EXPECT_NEAR(forwardDot, adjointDot, 2e-6);
        const auto zero = model.computeDose(std::vector<std::vector<float>>(3));
        for (double value : zero) EXPECT_DOUBLE_EQ(value, 0.0);
    };
    verify(gpu);
    auto streamed = cpu;
    const size_t budget = gpu.cudaMemoryBytes() - records.size() * sizeof(FluenceDoseModel::VoxelRecord);
    ASSERT_TRUE(streamed.enableCuda(budget, true));
    EXPECT_FALSE(streamed.cudaGeometryResident());
    verify(streamed);
    EXPECT_THROW(streamed.enableCuda(1, true), std::runtime_error);
    EXPECT_FALSE(streamed.usesCuda());
#endif
}

TEST(CudaFluenceTest, WaterCubeDoseAndBixelGradientMatchCpu) {
#ifndef OPTIRAD_HAS_CUDA
    GTEST_SKIP() << "CUDA backend not compiled";
#else
    ModelFixture fixture;
    FluenceModelOptions options;
    options.backend = DoseComputeBackend::CPU;
    auto cpu = std::make_shared<FluenceDoseModel>(fixture.engine.buildFluenceModel(
        *fixture.built.plan, *fixture.built.stf, *fixture.patient, *fixture.grid, options));
    options.backend = DoseComputeBackend::CUDA;
    auto gpu = std::make_shared<FluenceDoseModel>(fixture.engine.buildFluenceModel(
        *fixture.built.plan, *fixture.built.stf, *fixture.patient, *fixture.grid, options));
    ASSERT_TRUE(gpu->usesCuda());
    FluenceDoseOperator reference(cpu, fixture.built.stf), accelerated(gpu, fixture.built.stf);
    const auto cpuDose = reference.computeDose(fixture.weights);
    const auto gpuDose = accelerated.computeDose(fixture.weights);
    const double peak = *std::max_element(cpuDose.begin(), cpuDose.end());
    for (size_t voxel = 0; voxel < cpuDose.size(); ++voxel)
        EXPECT_NEAR(gpuDose[voxel], cpuDose[voxel], 2e-6 * peak);
    std::vector<double> voxelGradient(cpuDose.size());
    for (size_t voxel = 0; voxel < voxelGradient.size(); ++voxel)
        voxelGradient[voxel] = std::sin(0.013 * voxel);
    std::vector<double> cpuGradient(fixture.weights.size(), 0.0), gpuGradient = cpuGradient;
    reference.accumulateTransposeProduct(voxelGradient, cpuGradient);
    accelerated.accumulateTransposeProduct(voxelGradient, gpuGradient);
    for (size_t bixel = 0; bixel < cpuGradient.size(); ++bixel)
        EXPECT_NEAR(gpuGradient[bixel], cpuGradient[bixel], 2e-5 * std::max(1.0, std::abs(cpuGradient[bixel])));

    for (int type = 0; type < 5; ++type) {
        ObjectiveSpec objectiveSpec;
        objectiveSpec.structurePattern = "PTV";
        objectiveSpec.typeIdx = type;
        objectiveSpec.doseValue = peak * 0.4;
        objectiveSpec.volumePct = 35.0;
        objectiveSpec.weight = 3.0;
        auto built = ObjectiveBuilder::buildFromSpecs({objectiveSpec}, *fixture.patient, fixture.patient->getGrid(), *fixture.grid);
        ASSERT_EQ(built.ptrs.size(), 1u);
        DoseObjectiveSpec description;
        ASSERT_TRUE(built.ptrs.front()->gpuDescription(description));
        double actualValue = 0.0;
        std::vector<double> actualGradient;
        ASSERT_TRUE(accelerated.evaluateObjectives(fixture.weights, {description}, 0.0, 1.04, 0.0,
            actualValue, actualGradient));
        const double expectedValue = built.ptrs.front()->compute(cpuDose);
        EXPECT_NEAR(actualValue, expectedValue, 1e-5 * std::max(1.0, std::abs(expectedValue))) << "type " << type;
        const auto voxelDerivative = built.ptrs.front()->gradient(cpuDose);
        std::vector<double> expectedGradient(fixture.weights.size(), 0.0);
        reference.accumulateTransposeProduct(voxelDerivative, expectedGradient);
        for (size_t bixel = 0; bixel < expectedGradient.size(); ++bixel)
            EXPECT_NEAR(actualGradient[bixel], expectedGradient[bixel],
                1e-4 * std::max(1.0, std::abs(expectedGradient[bixel]))) << "type " << type;
    }

    for (bool consistent : {false, true}) {
        const double prescription = peak * 0.25, penalty = 2.0, threshold = 1.04;
        double expectedValue = 0.0, actualValue = 0.0;
        std::vector<double> derivative(cpuDose.size(), 0.0), expected(fixture.weights.size(), 0.0), actual;
        for (size_t voxel = 0; voxel < cpuDose.size(); ++voxel) {
            const double excess = cpuDose[voxel] - prescription * threshold;
            if (excess <= 0.0) continue;
            const double normalized = excess / prescription;
            expectedValue += penalty * normalized * normalized * (1.0 + 10.0 * normalized) * excess * excess;
            derivative[voxel] = consistent ? penalty * excess * normalized * normalized * (4.0 + 50.0 * normalized) :
                2.0 * penalty * normalized * (1.0 + 15.0 * normalized) * excess;
        }
        reference.accumulateTransposeProduct(derivative, expected);
        ASSERT_TRUE(accelerated.evaluateObjectives(fixture.weights, {}, prescription, threshold, penalty,
            actualValue, actual, consistent));
        EXPECT_NEAR(actualValue, expectedValue, 1e-5 * std::max(1.0, expectedValue));
        for (size_t bixel = 0; bixel < expected.size(); ++bixel)
            EXPECT_NEAR(actual[bixel], expected[bixel], 1e-4 * std::max(1.0, std::abs(expected[bixel])));
        ASSERT_TRUE(accelerated.evaluateObjectives(fixture.weights, {}, prescription, threshold, penalty,
            actualValue, actual, consistent, false));
        EXPECT_TRUE(actual.empty());
        EXPECT_NEAR(actualValue, expectedValue, 1e-5 * std::max(1.0, expectedValue));
    }

    ObjectiveSpec spec;
    spec.structurePattern = "PTV";
    spec.typeIdx = 0;
    spec.doseValue = 60.0;
    spec.weight = 1.0;
    auto objectives = ObjectiveBuilder::buildFromSpecs({spec}, *fixture.patient, fixture.patient->getGrid(), *fixture.grid);
    LBFGSOptimizer cpuOptimizer, gpuOptimizer;
    for (auto* optimizer : {&cpuOptimizer, &gpuOptimizer}) {
        optimizer->setMaxIterations(20);
        optimizer->setMaxFluence(1000.0);
        optimizer->setHotspotPenalty(0.0);
    }
    const auto cpuResult = cpuOptimizer.optimize(reference, objectives.ptrs, {});
    const auto gpuResult = gpuOptimizer.optimize(accelerated, objectives.ptrs, {});
    EXPECT_NEAR(gpuResult.finalObjective, cpuResult.finalObjective,
        1e-3 * std::max(1.0, std::abs(cpuResult.finalObjective)));

    auto milliseconds = [&](const IDoseOperator& doseOperator) {
        const auto start = std::chrono::steady_clock::now();
        for (int repeat = 0; repeat < 5; ++repeat) {
            (void)doseOperator.computeDose(fixture.weights);
            std::vector<double> gradient(fixture.weights.size(), 0.0);
            doseOperator.accumulateTransposeProduct(voxelGradient, gradient);
        }
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count() / 5.0;
    };
    std::cout << "Forward+adjoint: CPU " << milliseconds(reference) << " ms, CUDA "
              << milliseconds(accelerated) << " ms (3-beam water cube)\n";
#endif
}

TEST(CudaRayTracerTest, MatchesCpuOnNonuniformDensityAndNoncubicGrid) {
#ifndef OPTIRAD_HAS_CUDA
    GTEST_SKIP() << "CUDA backend not compiled";
#else
    Grid grid;
    grid.setDimensions(19, 23, 17);
    grid.setSpacing(2.0, 3.0, 4.0);
    grid.setOrigin({-20.0, -30.0, -32.0});
    std::vector<double> density(grid.getNumVoxels());
    for (size_t voxel = 0; voxel < density.size(); ++voxel) density[voxel] = 0.2 + (voxel % 23) * 0.07;
    CudaRayTracer tracer(grid, density.data());
    std::vector<Vec3> targets;
    for (int index = 0; index < 200; ++index)
        targets.push_back({-35.0 + (index % 17) * 4.3, -33.0 + (index % 19) * 3.7, -38.0 + (index % 13) * 5.1});
    for (const Vec3 source : {Vec3{0.0, -1000.0, 0.0}, Vec3{500.0, 700.0, -600.0}, Vec3{0.0, 0.0, 0.0}}) {
        const auto actual = tracer.depths(source, targets);
        for (size_t ray = 0; ray < targets.size(); ++ray) {
            const auto cpu = SiddonRayTracer::traceRadDepth(source, targets[ray], grid, density.data());
            EXPECT_NEAR(actual[ray], cpu.empty() ? 0.0 : cpu.back().second, 1e-7) << "ray " << ray;
        }
    }
#endif
}

TEST(CudaDirectDoseTest, MatchesCpuReferenceDeposition) {
#ifndef OPTIRAD_HAS_CUDA
    GTEST_SKIP() << "CUDA backend not compiled";
#else
    ModelFixture fixture;
    const auto begin = std::chrono::steady_clock::now();
    const auto cpu = fixture.engine.calculateDoseDirect(*fixture.built.plan, *fixture.built.stf,
        *fixture.patient, *fixture.grid, fixture.weights);
    const auto middle = std::chrono::steady_clock::now();
    auto options = fixture.engine.getOptions();
    options.gpuRayTracing = true;
    options.requireGpuRayTracing = true;
    fixture.engine.setOptions(options);
    const auto gpu = fixture.engine.calculateDoseDirect(*fixture.built.plan, *fixture.built.stf,
        *fixture.patient, *fixture.grid, fixture.weights);
    const auto end = std::chrono::steady_clock::now();
    ASSERT_EQ(cpu.size(), gpu.size());
    const double tolerance = 1e-7 * std::max(1.0, cpu.getMax());
    for (size_t voxel = 0; voxel < cpu.size(); ++voxel)
        EXPECT_NEAR(gpu.data()[voxel], cpu.data()[voxel], tolerance);
    std::cout << "Direct dose CPU " << std::chrono::duration<double>(middle - begin).count()
              << " s, CUDA " << std::chrono::duration<double>(end - middle).count() << " s\n";
#endif
}

TEST(ApertureDoseOperatorTest, MatchesExpandedDoseTransposeAndCudaObjectives) {
    ModelFixture fixture;
    FluenceModelOptions options;
    options.backend = DoseComputeBackend::CPU;
    auto model = std::make_shared<FluenceDoseModel>(fixture.engine.buildFluenceModel(
        *fixture.built.plan, *fixture.built.stf, *fixture.patient, *fixture.grid, options));
    auto bixels = std::make_shared<FluenceDoseOperator>(model, fixture.built.stf);
    ApertureDoseOperator apertures(bixels, *fixture.built.stf, fixture.weights);
    const std::vector<double> weights{0.7, 0.0, 1.4};
    auto expanded = fixture.weights;
    for (size_t beam = 0; beam < weights.size(); ++beam)
        for (size_t bixel = fixture.offsets[beam]; bixel < fixture.offsets[beam + 1]; ++bixel)
            expanded[bixel] *= weights[beam];
    const auto expected = bixels->computeDose(expanded);
    const auto actual = apertures.computeDose(weights);
    ASSERT_EQ(expected.size(), actual.size());
    for (size_t voxel = 0; voxel < actual.size(); ++voxel) EXPECT_DOUBLE_EQ(actual[voxel], expected[voxel]);
    std::vector<double> voxelGradient(actual.size());
    for (size_t voxel = 0; voxel < actual.size(); ++voxel) voxelGradient[voxel] = std::sin(0.013 * voxel);
    std::vector<double> gradient(weights.size(), 0.0);
    apertures.accumulateTransposeProduct(voxelGradient, gradient);
    double lhs = 0.0, rhs = 0.0;
    for (size_t voxel = 0; voxel < actual.size(); ++voxel) lhs += actual[voxel] * voxelGradient[voxel];
    for (size_t aperture = 0; aperture < weights.size(); ++aperture) rhs += weights[aperture] * gradient[aperture];
    EXPECT_NEAR(lhs, rhs, 1e-4 * std::max(1.0, std::abs(lhs)));
    auto accumulated = gradient;
    apertures.accumulateTransposeProduct(voxelGradient, accumulated);
    for (size_t aperture = 0; aperture < weights.size(); ++aperture)
        EXPECT_NEAR(accumulated[aperture], 2.0 * gradient[aperture], 1e-8);
    EXPECT_THROW(apertures.computeDose({1.0}), std::invalid_argument);

#ifdef OPTIRAD_HAS_CUDA
    auto gpuModel = std::make_shared<FluenceDoseModel>(*model);
    ASSERT_TRUE(gpuModel->enableCuda(0, true));
    auto gpuBixels = std::make_shared<FluenceDoseOperator>(gpuModel, fixture.built.stf);
    ApertureDoseOperator gpu(gpuBixels, *fixture.built.stf, fixture.weights);
    DoseObjectiveSpec objective;
    objective.type = 0;
    objective.threshold = 2.0;
    objective.weight = 3.0;
    for (size_t voxel = 0; voxel < actual.size(); voxel += 17) objective.indices.push_back(voxel);
    double expectedValue = 0.0, gpuValue = 0.0;
    std::vector<double> derivative(actual.size(), 0.0), expectedGradient(weights.size(), 0.0), gpuGradient;
    for (size_t voxel : objective.indices) {
        const double difference = actual[voxel] - objective.threshold;
        expectedValue += objective.weight * difference * difference / objective.indices.size();
        derivative[voxel] = 2.0 * objective.weight * difference / objective.indices.size();
    }
    apertures.accumulateTransposeProduct(derivative, expectedGradient);
    ASSERT_TRUE(gpu.evaluateObjectives(weights, {objective}, 0.0, 1.04, 0.0, gpuValue, gpuGradient));
    EXPECT_NEAR(gpuValue, expectedValue, 1e-5 * std::max(1.0, std::abs(expectedValue)));
    ASSERT_EQ(gpuGradient.size(), weights.size());
    for (size_t aperture = 0; aperture < weights.size(); ++aperture)
        EXPECT_NEAR(gpuGradient[aperture], expectedGradient[aperture], 1e-5 * std::max(1.0, std::abs(expectedGradient[aperture])));
    ASSERT_TRUE(gpu.evaluateObjectives(weights, {objective}, 0.0, 1.04, 0.0, gpuValue, gpuGradient, false, false));
    EXPECT_TRUE(gpuGradient.empty());
#endif
}

} // namespace optirad::tests
