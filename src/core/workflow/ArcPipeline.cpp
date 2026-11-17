#include "ArcPipeline.hpp"

#include "dose/DoseEngineFactory.hpp"
#include "dose/FluenceDoseOperator.hpp"
#include "dose/engines/PencilBeamEngine.hpp"
#include "optimization/OptimizationLogging.hpp"
#include "sequencing/ApertureFluence.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "utils/Logger.hpp"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <chrono>
#include <iomanip>
#include <sstream>

namespace optirad {

namespace {

// Aperture of the adjacent control point of the same arc, or null at the ends of an arc.
const Aperture* neighbour(const ArcSequenceResult& sequence, size_t index, int step) {
    if ((step < 0 && index == 0) || (step > 0 && index + 1 >= sequence.controlPoints.size())) return nullptr;
    const auto& other = sequence.controlPoints[step < 0 ? index - 1 : index + 1];
    const auto& cp = sequence.controlPoints[index];
    if (other.geometry.arcIndex != cp.geometry.arcIndex ||
        other.aperture.bankA.size() != cp.aperture.bankA.size()) return nullptr;
    return &other.aperture;
}

// One run of the leaf/MU optimizer on `result.sequence`. The progressive optimizer calls it once per resolution.
struct RefineOptions {
    int stage = 1;
    int stages = 1;
    ArcProgressCallback progressCallback;
    int iterations = 0;
    double maxMoveMm = 0.0;         // largest leaf tip move per iteration; 0 = 4 x leafStepMm
    double* muLambda = nullptr;     // objective units per MU: < 0 on entry = derive from the initial objective
    double* smoothingReference = nullptr;   // objective scale of the MU smoothing: < 0 on entry = initial objective
    double muPenalty = 0.0;         // fraction of the starting objective that the total MU costs
    bool reportInitial = true;
    bool finalStage = true;         // snap to the MLC resolution and report timing at the end
};

void refineLeaves(ArcPipelineResult& result, const Plan& plan, const PatientData& patient,
                  const Grid& grid, const BuiltObjectives& objectives, const ArcOptimizationConfig& config,
                  std::shared_ptr<FluenceDoseModel> preparedModel, const RefineOptions& refine) {
    const auto startTime = std::chrono::steady_clock::now();
    const int initialAcceptedSteps = result.leafIterations;
    if (refine.iterations < 1 || config.leafStepMm <= 0.0 || !config.engineOptions.divergentBixels)
        throw std::invalid_argument("Leaf optimization requires positive iterations/step and divergent bixels");
    PencilBeamEngine engine;
    engine.setOptions(config.engineOptions);
    FluenceModelOptions options;
    options.pixelMm = config.fluencePixelMm;
    options.backend = config.computeBackend;
    options.gpuMemoryBudgetMiB = config.gpuMemoryBudgetMiB;
    auto model = std::move(preparedModel);
    if (!model) model = std::make_shared<FluenceDoseModel>(engine.buildFluenceModel(
        plan, *result.controlPointStf, patient, grid, options));
    FluenceDoseOperator doseOperator(model, result.controlPointStf);
    Logger::info(std::string("DAO leaf refinement: dose/objective backend=") +
                 (model->usesCuda() ? "CUDA" : "CPU") +
                 ", control points=" + std::to_string(result.sequence.controlPoints.size()));
    const MachineConstraints limits = ArcPipeline::deliveryConstraints(plan);
    const auto& mlc = plan.getMachine().getGeometry();
    const double resolution = mlc.leafPositionResolution;
    auto snap = [resolution](double position) {
        return resolution > 0.0 ? std::round(position / resolution) * resolution : position;
    };
    for (auto& cp : result.sequence.controlPoints) {
        auto& aperture = cp.aperture;
        const double low = mlc.maxLeafTravel > 0.0 ? std::max(-mlc.jawX1Max, -mlc.maxLeafTravel) : -mlc.jawX1Max;
        const double high = mlc.maxLeafTravel > 0.0 ? std::min(mlc.jawX2Max, mlc.maxLeafTravel) : mlc.jawX2Max;
        for (double& tip : aperture.bankA) tip = snap(std::clamp(tip, low, high));
        for (double& tip : aperture.bankB) tip = snap(std::clamp(tip, low, high));
        if (!mlc.interdigitation && !aperture.bankA.empty()) {
            double maxLeft = low, minRight = high;
            for (size_t pair = 0; pair < aperture.bankA.size(); ++pair) {
                if (aperture.bankA[pair] >= aperture.bankB[pair]) continue;
                maxLeft = std::max(maxLeft, aperture.bankA[pair]);
                minRight = std::min(minRight, aperture.bankB[pair]);
            }
            const double meeting = snap(maxLeft <= minRight ? std::clamp(0.0, maxLeft, minRight) :
                0.5 * (maxLeft + minRight));
            for (size_t pair = 0; pair < aperture.bankA.size(); ++pair) {
                if (aperture.bankA[pair] >= aperture.bankB[pair]) {
                    aperture.bankA[pair] = aperture.bankB[pair] = meeting;
                } else if (maxLeft > minRight) {
                    aperture.bankA[pair] = std::min(aperture.bankA[pair], meeting);
                    aperture.bankB[pair] = std::max(aperture.bankB[pair], meeting);
                }
            }
        }
    }
    ArcSequencer::computeTiming(result.sequence, limits);
    for (auto& cp : result.sequence.controlPoints) {
        if (cp.mu > 0.0) cp.mu = std::max(cp.mu, limits.monitorUnitRate[0] * cp.timeSec);
        cp.aperture.weight = cp.mu;
    }
    ArcSequencer::computeTiming(result.sequence, limits);
    const auto seed = result.sequence;
    // The gantry may slow down for a control point that needs more MU than the maximum dose rate delivers at full
    // speed (as on the machine and in Eclipse); the leaf travel limits keep the full-speed time (conservative).
    const double slowestGantry = limits.gantryRotationSpeed[0] > 0.0 ? limits.gantryRotationSpeed[0] :
        0.1 * (limits.gantryRotationSpeed[1] > 0.0 ? limits.gantryRotationSpeed[1] : 6.0);
    auto maxMuOf = [&](const DeliveredControlPoint& cp) {
        return limits.monitorUnitRate[1] * std::max(cp.timeSec, cp.geometry.width / slowestGantry);
    };
    const size_t count = seed.controlPoints.size();

    // With `quantize` leaf tips are snapped to the MLC position resolution; the iterations run on continuous positions.
    auto project = [&](ArcSequenceResult& sequence, bool quantize) {
        auto place = [&](double position) { return quantize ? snap(position) : position; };
        const double step = quantize ? resolution : 0.0;
        for (size_t index = 0; index < count; ++index) {
            auto& cp = sequence.controlPoints[index];
            const auto& original = seed.controlPoints[index];
            const double duration = original.timeSec;
            double low = -mlc.jawX1Max, high = mlc.jawX2Max;
            if (mlc.maxLeafTravel > 0.0) {
                low = std::max(low, -mlc.maxLeafTravel);
                high = std::min(high, mlc.maxLeafTravel);
            }
            auto& aperture = cp.aperture;
            double minMeeting = low, maxMeeting = high;
            for (size_t pair = 0; pair < aperture.bankA.size(); ++pair) {
                double lowA = low, highA = high, lowB = low, highB = high;
                if (index > 0 && cp.geometry.arcIndex == sequence.controlPoints[index - 1].geometry.arcIndex) {
                    const auto& previous = sequence.controlPoints[index - 1];
                    const double travel = duration * limits.leafSpeed[1];
                    lowA = std::max(lowA, previous.aperture.bankA[pair] - travel);
                    highA = std::min(highA, previous.aperture.bankA[pair] + travel);
                    lowB = std::max(lowB, previous.aperture.bankB[pair] - travel);
                    highB = std::min(highB, previous.aperture.bankB[pair] + travel);
                }
                if (step > 0.0) {
                    lowA = std::ceil((lowA - 1e-9) / step) * step;
                    highA = std::floor((highA + 1e-9) / step) * step;
                    lowB = std::ceil((lowB - 1e-9) / step) * step;
                    highB = std::floor((highB + 1e-9) / step) * step;
                }
                if (lowA > highA || lowB > highB) return false;
                aperture.bankA[pair] = std::clamp(place(aperture.bankA[pair]), lowA, highA);
                aperture.bankB[pair] = std::clamp(place(aperture.bankB[pair]), lowB, highB);
                minMeeting = std::max(minMeeting, lowA);
                maxMeeting = std::min(maxMeeting, highB);
                if (aperture.bankA[pair] > aperture.bankB[pair]) {
                    if (lowA > highB) return false;
                    const double meeting = std::clamp(place(0.5 * (aperture.bankA[pair] + aperture.bankB[pair])),
                        lowA, highB);
                    aperture.bankA[pair] = std::min(aperture.bankA[pair], meeting);
                    aperture.bankB[pair] = std::max(aperture.bankB[pair], meeting);
                }
            }
            if (!mlc.interdigitation && !aperture.bankA.empty()) {
                if (minMeeting > maxMeeting) return false;
                const double leftmost = *std::max_element(aperture.bankA.begin(), aperture.bankA.end());
                const double rightmost = *std::min_element(aperture.bankB.begin(), aperture.bankB.end());
                if (leftmost > rightmost) {
                    const double meeting = std::clamp(place(0.5 * (leftmost + rightmost)), minMeeting, maxMeeting);
                    for (double& tip : aperture.bankA) tip = std::min(tip, meeting);
                    for (double& tip : aperture.bankB) tip = std::max(tip, meeting);
                }
            }
            if (duration <= 0.0 || limits.leafSpeed[1] <= 0.0 ||
                cp.geometry.width / duration < limits.gantryRotationSpeed[0]) return false;
            const double minMU = original.mu > 0.0 ? limits.monitorUnitRate[0] * duration : 0.0;
            const double maxMU = maxMuOf(original);
            if (minMU > maxMU) return false;
            cp.mu = std::clamp(cp.mu, minMU, maxMU);
            aperture.weight = cp.mu;
            cp.timeSec = duration;
            cp.gantrySpeed = cp.geometry.width / duration;
            cp.doseRate = cp.mu / duration;
        }
        return true;
    };

    std::vector<DoseObjectiveSpec> descriptions;
    bool supported = true;
    for (const auto* objective : objectives.ptrs) {
        DoseObjectiveSpec description;
        if (!objective->gpuDescription(description)) { supported = false; break; }
        descriptions.push_back(std::move(description));
    }
    double muSmooth = 0.0;   // objective units per MU^2, set from config.muSmoothing once the MU scale is known
    auto evaluate = [&](const ArcSequenceResult& sequence, std::vector<double>* gradient) {
        const auto weights = ArcPipeline::bixelWeights(*result.controlPointStf, sequence, config.integrateLeafMotion);
        double value = 0.0;
        std::vector<double> unusedGradient;
        if (!supported || !doseOperator.evaluateObjectives(weights, descriptions, value,
            gradient ? *gradient : unusedGradient, gradient != nullptr)) {
        const auto dose = doseOperator.computeDose(weights);
        std::vector<double> voxelGradient(gradient ? dose.size() : 0, 0.0);
        for (const auto* objective : objectives.ptrs) {
            value += objective->compute(dose);
            if (gradient) {
                const auto contribution = objective->gradient(dose);
                for (size_t voxel = 0; voxel < dose.size(); ++voxel) voxelGradient[voxel] += contribution[voxel];
            }
        }
        if (gradient) {
            gradient->assign(weights.size(), 0.0);
            doseOperator.accumulateTransposeProduct(voxelGradient, *gradient);
        }
        }
        for (size_t index = 1; index < count; ++index) {
            const auto& cp = sequence.controlPoints[index];
            const auto& previous = sequence.controlPoints[index - 1];
            if (cp.geometry.arcIndex != previous.geometry.arcIndex) continue;
            const double difference = cp.mu - previous.mu;
            value += muSmooth * difference * difference;
        }
        return value;
    };

    if (!project(result.sequence, false))
        throw std::runtime_error("VMAT seed cannot satisfy the configured leaf/MU bounds in its time slots");

    // MU are optimized in units of 10 % of their seed value, tips in mm.
    std::vector<double> muScale(count);
    double seedMu = 0.0;
    for (size_t index = 0; index < count; ++index) {
        muScale[index] = 0.1 * std::max(seed.controlPoints[index].mu, 1e-6);
        seedMu += seed.controlPoints[index].mu;
    }
    double ownLambda = 0.0;
    double& muLambda = refine.muLambda ? *refine.muLambda : ownLambda;   // objective units per MU
    const bool deriveLambda = muLambda < 0.0;
    if (deriveLambda) muLambda = 0.0;

    // Componentwise bounds: jaw/travel limits of the tips and the dose-rate limits of the MU of each control point.
    double tipLow = -mlc.jawX1Max, tipHigh = mlc.jawX2Max;
    if (mlc.maxLeafTravel > 0.0) {
        tipLow = std::max(tipLow, -mlc.maxLeafTravel);
        tipHigh = std::min(tipHigh, mlc.maxLeafTravel);
    }
    std::vector<double> muLow(count, 0.0), muHigh(count), travelLimit(count, 0.0);
    for (size_t index = 0; index < count; ++index) {
        const auto& original = seed.controlPoints[index];
        muLow[index] = original.mu > 0.0 ? limits.monitorUnitRate[0] * original.timeSec : 0.0;
        muHigh[index] = maxMuOf(original);
        travelLimit[index] = limits.leafSpeed[1] * original.timeSec;
    }
    auto boxClamp = [&](ArcSequenceResult& sequence) {
        for (size_t index = 0; index < count; ++index) {
            auto& cp = sequence.controlPoints[index];
            for (double& tip : cp.aperture.bankA) tip = std::clamp(tip, tipLow, tipHigh);
            for (double& tip : cp.aperture.bankB) tip = std::clamp(tip, tipLow, tipHigh);
            cp.mu = std::clamp(cp.mu, muLow[index], std::max(muLow[index], muHigh[index]));
        }
    };
    // The leaf speed limit and the order of the leaves of a pair are quadratic penalties (a projection that
    // couples the control points turns descent directions into ascent directions); `project` repairs the result.
    double penaltyWeight = 0.0;

    // Objective (with the MU penalty) and, optionally, its gradient in the scaled variables, stored in a copy of
    // the sequence: tip entries are mm gradients, `mu` is the gradient with respect to MU / muScale.
    auto objective = [&](const ArcSequenceResult& sequence, ArcSequenceResult* direction) {
        std::vector<double> gradient;
        double value = evaluate(sequence, direction ? &gradient : nullptr);
        for (const auto& cp : sequence.controlPoints) value += muLambda * cp.mu;
        if (!direction) return value;

        *direction = sequence;
        size_t offset = 0;
        for (auto& delta : direction->controlPoints) {
            std::fill(delta.aperture.bankA.begin(), delta.aperture.bankA.end(), 0.0);
            std::fill(delta.aperture.bankB.begin(), delta.aperture.bankB.end(), 0.0);
            delta.mu = 0.0;
        }
        for (size_t index = 0; index < count; ++index) {
            const auto& beam = *result.controlPointStf->getBeam(index);
            const auto& cp = sequence.controlPoints[index];
            auto& delta = direction->controlPoints[index];
            const auto* before = neighbour(sequence, index, -1);
            const auto* after = neighbour(sequence, index, 1);
            const auto open = config.integrateLeafMotion ?
                ApertureFluence::integratedOpenFractions(beam, before, cp.aperture, after, cp.leafPairBoundariesZ) :
                ApertureFluence::bixelOpenFractions(beam, cp.aperture, cp.leafPairBoundariesZ);
            const double width = beam.getBixelWidth();
            std::vector<double> rayWeight(beam.getNumOfRays());
            for (size_t ray = 0; ray < beam.getNumOfRays(); ++ray) {
                delta.mu += gradient[offset + ray] * open[ray];
                rayWeight[ray] = cp.mu * gradient[offset + ray];
            }
            if (config.integrateLeafMotion) {
                ApertureFluence::accumulateIntegratedTipGradient(beam, before, cp.aperture, after,
                    cp.leafPairBoundariesZ, rayWeight, before ? &direction->controlPoints[index - 1].aperture : nullptr,
                    delta.aperture, after ? &direction->controlPoints[index + 1].aperture : nullptr);
            } else {
                for (size_t ray = 0; ray < beam.getNumOfRays(); ++ray) {
                    const auto position = beam.getRay(ray)->getRayPosBev();
                    const double left = position[0] - width * 0.5, right = left + width;
                    const double bottom = position[2] - width * 0.5, top = bottom + width;
                    for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
                        const double overlap = std::max(0.0, std::min(top, cp.leafPairBoundariesZ[pair + 1]) -
                            std::max(bottom, cp.leafPairBoundariesZ[pair]));
                        const double scale = rayWeight[ray] * overlap / (width * width);
                        const double bankA = cp.aperture.bankA[pair], bankB = cp.aperture.bankB[pair];
                        if (bankA >= left && bankA < right && bankA <= bankB) delta.aperture.bankA[pair] -= scale;
                        if (bankB > left && bankB <= right && bankB >= bankA) delta.aperture.bankB[pair] += scale;
                    }
                }
            }
            offset += beam.getTotalNumOfBixels();
        }
        for (size_t index = 0; index < count; ++index) {
            auto& delta = direction->controlPoints[index];
            delta.mu += muLambda;
            if (index == 0) continue;
            const auto& cp = sequence.controlPoints[index];
            const auto& previous = sequence.controlPoints[index - 1];
            if (cp.geometry.arcIndex != previous.geometry.arcIndex) continue;
            const double derivative = 2.0 * muSmooth * (cp.mu - previous.mu);
            delta.mu += derivative;
            direction->controlPoints[index - 1].mu -= derivative;
        }
        if (penaltyWeight > 0.0) {
            for (size_t index = 0; index < count; ++index) {
                const auto& cp = sequence.controlPoints[index];
                auto& delta = direction->controlPoints[index].aperture;
                for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
                    const double order = cp.aperture.bankA[pair] - cp.aperture.bankB[pair];
                    if (order > 0.0) {
                        value += penaltyWeight * order * order;
                        delta.bankA[pair] += 2.0 * penaltyWeight * order;
                        delta.bankB[pair] -= 2.0 * penaltyWeight * order;
                    }
                }
                if (index == 0 || cp.geometry.arcIndex != sequence.controlPoints[index - 1].geometry.arcIndex) continue;
                const auto& previous = sequence.controlPoints[index - 1];
                auto& deltaPrevious = direction->controlPoints[index - 1].aperture;
                for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
                    const double moveA = cp.aperture.bankA[pair] - previous.aperture.bankA[pair];
                    const double moveB = cp.aperture.bankB[pair] - previous.aperture.bankB[pair];
                    for (int bank = 0; bank < 2; ++bank) {
                        const double move = bank == 0 ? moveA : moveB;
                        const double excess = std::abs(move) - travelLimit[index];
                        if (excess <= 0.0) continue;
                        value += penaltyWeight * excess * excess;
                        const double slope = 2.0 * penaltyWeight * excess * (move > 0.0 ? 1.0 : -1.0);
                        (bank == 0 ? delta.bankA : delta.bankB)[pair] += slope;
                        (bank == 0 ? deltaPrevious.bankA : deltaPrevious.bankB)[pair] -= slope;
                    }
                }
            }
        }
        for (size_t index = 0; index < count; ++index)
            direction->controlPoints[index].mu *= muScale[index];
        return value;
    };

    result.initialLeafObjective = evaluate(result.sequence, nullptr);
    if (deriveLambda && refine.muPenalty > 0.0 && seedMu > 0.0)
        muLambda = refine.muPenalty * result.initialLeafObjective / seedMu;
    if (refine.smoothingReference && config.muSmoothing > 0.0 && seedMu > 0.0) {
        double& reference = *refine.smoothingReference;
        if (reference < 0.0) reference = result.initialLeafObjective;
        size_t neighbours = 0;
        for (size_t index = 1; index < count; ++index)
            if (seed.controlPoints[index].geometry.arcIndex == seed.controlPoints[index - 1].geometry.arcIndex)
                ++neighbours;
        const double meanMu = seedMu / static_cast<double>(count);
        if (neighbours > 0) muSmooth = config.muSmoothing * reference / (static_cast<double>(neighbours) * meanMu * meanMu);
    }

    // Projected L-BFGS over all leaf tips (mm) and MU (in units of muScale). The speed and MU-rate bounds couple
    // neighbouring control points, so they are enforced by `project`; a step that is not a descent direction after
    // projection falls back to a (scaled) gradient step with a cleared history.
    const size_t pairsPerCp = result.sequence.controlPoints.empty() ? 0 : result.sequence.controlPoints.front().aperture.bankA.size();
    const size_t perCp = 2 * pairsPerCp + 1;
    const size_t dim = perCp * count;
    auto pack = [&](const ArcSequenceResult& sequence, bool fromGradient) {
        std::vector<double> x(dim);
        for (size_t index = 0; index < count; ++index) {
            const auto& cp = sequence.controlPoints[index];
            double* out = &x[index * perCp];
            for (size_t pair = 0; pair < pairsPerCp; ++pair) {
                out[pair] = cp.aperture.bankA[pair];
                out[pairsPerCp + pair] = cp.aperture.bankB[pair];
            }
            out[2 * pairsPerCp] = fromGradient ? cp.mu : cp.mu / muScale[index];
        }
        return x;
    };
    auto unpack = [&](const std::vector<double>& x, ArcSequenceResult& sequence) {
        for (size_t index = 0; index < count; ++index) {
            auto& cp = sequence.controlPoints[index];
            const double* in = &x[index * perCp];
            for (size_t pair = 0; pair < pairsPerCp; ++pair) {
                cp.aperture.bankA[pair] = in[pair];
                cp.aperture.bankB[pair] = in[pairsPerCp + pair];
            }
            cp.mu = in[2 * pairsPerCp] * muScale[index];
        }
    };
    auto dot = [](const std::vector<double>& a, const std::vector<double>& b) {
        double sum = 0.0;
        for (size_t i = 0; i < a.size(); ++i) sum += a[i] * b[i];
        return sum;
    };

    ArcSequenceResult current = result.sequence, direction, trial, trialDirection;
    double value = objective(current, &direction);
    {   // 1 mm of speed violation costs as much as ten times the steepest tip gradient
        double steepest = 0.0;
        for (const auto& delta : direction.controlPoints) {
            for (double v : delta.aperture.bankA) steepest = std::max(steepest, std::abs(v));
            for (double v : delta.aperture.bankB) steepest = std::max(steepest, std::abs(v));
        }
        penaltyWeight = 10.0 * steepest;
        value = objective(current, &direction);
    }
    ArcSequenceResult best = current;
    double bestValue = value;
    std::vector<double> bestHistory{value};
    std::vector<double> x = pack(current, false), g = pack(direction, true);
    std::deque<std::vector<double>> sHistory, yHistory;
    const size_t memory = 10;
    const double maxMove = std::max(refine.maxMoveMm > 0.0 ? refine.maxMoveMm : 4.0 * config.leafStepMm, 1e-3);
    const double maxMuMove = 4.0;   // in muScale units (40 % of the seed MU)
    {
        std::ostringstream out;
        out << std::setprecision(5) << "DAO settings: iterations=" << refine.iterations
            << ", leaf step=" << config.leafStepMm << " mm, max move=" << maxMove << " mm"
            << "\n  Penalties: MU smooth=" << muSmooth << " (relative " << config.muSmoothing << ")"
            << ", MU total=" << muLambda << " (seed=" << refine.muPenalty << "), leaf=" << penaltyWeight
            << "\n  Bounds: tips=[" << tipLow << ", " << tipHigh
            << "] mm, leaf<=" << limits.leafSpeed[1]
            << " mm/s, gantry<=" << limits.gantryRotationSpeed[1]
            << " deg/s\n  MU rate=[" << limits.monitorUnitRate[0] << ", " << limits.monitorUnitRate[1]
            << "] /s (1 MU = " << ArcPipeline::machineMuPerWeight(plan) << " machine MU/fraction), snap="
            << (refine.finalStage ? resolution : 0.0)
            << " mm, interdigitation=" << (mlc.interdigitation ? "allowed" : "forbidden");
        Logger::info(out.str());
    }
    auto report = [&](const std::string& context, const ArcSequenceResult& sequence,
                      DaoProgressInfo info, bool terminal, bool finished) {
        double mu = 0.0, smoothing = 0.0, orderPenalty = 0.0, speedPenalty = 0.0;
        for (size_t index = 0; index < count; ++index) {
            const auto& cp = sequence.controlPoints[index];
            mu += cp.mu;
            for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
                const double excess = std::max(0.0, cp.aperture.bankA[pair] - cp.aperture.bankB[pair]);
                orderPenalty += penaltyWeight * excess * excess;
            }
            if (index == 0 || cp.geometry.arcIndex != sequence.controlPoints[index - 1].geometry.arcIndex) continue;
            const auto& previous = sequence.controlPoints[index - 1];
            const double difference = cp.mu - previous.mu;
            smoothing += muSmooth * difference * difference;
            for (size_t pair = 0; pair < cp.aperture.bankA.size(); ++pair) {
                const double excessA = std::max(0.0,
                    std::abs(cp.aperture.bankA[pair] - previous.aperture.bankA[pair]) - travelLimit[index]);
                const double excessB = std::max(0.0,
                    std::abs(cp.aperture.bankB[pair] - previous.aperture.bankB[pair]) - travelLimit[index]);
                speedPenalty += penaltyWeight * (excessA * excessA + excessB * excessB);
            }
        }
        auto breakdown = computeOptimizationBreakdown(objectives.ptrs, doseOperator.computeDose(
            ArcPipeline::bixelWeights(*result.controlPointStf, sequence, config.integrateLeafMotion)),
            {{"MU smoothing", smoothing}, {"Total MU", muLambda * mu},
             {"Leaf order", orderPenalty}, {"Leaf speed", speedPenalty}});
        if (terminal) logOptimizationBreakdown(context, breakdown);
        if (refine.progressCallback) {
            info.totalMu = mu;
            info.objective = breakdown.total;
            info.elapsedSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
            ArcProgressSnapshot snapshot{refine.stage, refine.stages, refine.iterations, count,
                                         finished, info, std::move(breakdown)};
            result.progressHistory.push_back(snapshot);
            refine.progressCallback(snapshot);
        }
    };
    if (refine.reportInitial || refine.progressCallback)
        report("DAO initial", current, {}, refine.reportInitial, false);
    std::string stopReason = "iteration budget reached";
    int attemptedIterations = 0;
    int directionRestarts = 0, lineSearchRetries = 0;
    for (int iteration = 0; iteration < refine.iterations; ++iteration) {
        attemptedIterations = iteration + 1;
        double maxGradient = 0.0;
        for (double v : g) maxGradient = std::max(maxGradient, std::abs(v));
        if (maxGradient <= 0.0) {
            stopReason = "zero gradient";
            break;
        }

        // Variables on a bound that the step would leave stay where they are.
        auto freeze = [&](std::vector<double>& dir) {
            for (size_t k = 0; k < dim; ++k) {
                const size_t index = k / perCp, slot = k % perCp;
                const bool isMu = slot == perCp - 1;
                const double low = isMu ? muLow[index] / muScale[index] : tipLow;
                const double high = isMu ? std::max(muLow[index], muHigh[index]) / muScale[index] : tipHigh;
                if ((x[k] <= low + 1e-12 && dir[k] < 0.0) || (x[k] >= high - 1e-12 && dir[k] > 0.0)) dir[k] = 0.0;
            }
        };
        // The quasi-Newton step works on the free variables only (projected gradient); otherwise the curvature of
        // the active bounds (e.g. MU at the minimum dose rate) turns the direction uphill and forces restarts.
        std::vector<double> d(g);
        for (double& v : d) v = -v;
        freeze(d);
        for (double& v : d) v = -v;
        std::vector<double> alphas(sHistory.size());
        for (size_t i = sHistory.size(); i-- > 0;) {
            alphas[i] = dot(sHistory[i], d) / dot(yHistory[i], sHistory[i]);
            for (size_t k = 0; k < dim; ++k) d[k] -= alphas[i] * yHistory[i][k];
        }
        if (!sHistory.empty()) {
            const double gamma = dot(sHistory.back(), yHistory.back()) / dot(yHistory.back(), yHistory.back());
            for (double& v : d) v *= gamma;
        } else {
            for (double& v : d) v *= config.leafStepMm / maxGradient;
        }
        for (size_t i = 0; i < sHistory.size(); ++i) {
            const double beta = dot(yHistory[i], d) / dot(yHistory[i], sHistory[i]);
            for (size_t k = 0; k < dim; ++k) d[k] += (alphas[i] - beta) * sHistory[i][k];
        }
        for (double& v : d) v = -v;
        freeze(d);
        if (dot(d, g) >= 0.0) {   // not a descent direction: restart from the gradient
            sHistory.clear();
            yHistory.clear();
            for (size_t k = 0; k < dim; ++k) d[k] = -g[k] * config.leafStepMm / maxGradient;
            freeze(d);
            if (dot(d, g) >= 0.0) {
                stopReason = "stationary on the feasible box";
                break;
            }
            ++directionRestarts;
        }
        double biggestTip = 0.0, biggestMu = 0.0;
        for (size_t k = 0; k < dim; ++k) {
            double& biggest = (k % perCp == perCp - 1) ? biggestMu : biggestTip;
            biggest = std::max(biggest, std::abs(d[k]));
        }
        const double shrink = std::min({1.0, biggestTip > maxMove ? maxMove / biggestTip : 1.0,
                                        biggestMu > maxMuMove ? maxMuMove / biggestMu : 1.0});
        if (shrink < 1.0) for (double& v : d) v *= shrink;

        bool accepted = false;
        double step = 1.0;
        int lineSearchAttempts = 0;
        std::vector<double> xTrial(dim);
        for (int attempt = 0; attempt < 20; ++attempt, step *= 0.5) {
            ++lineSearchAttempts;
            for (size_t k = 0; k < dim; ++k) xTrial[k] = x[k] + step * d[k];
            trial = current;
            unpack(xTrial, trial);
            boxClamp(trial);
            xTrial = pack(trial, false);
            double slope = 0.0;
            for (size_t k = 0; k < dim; ++k) slope += g[k] * (xTrial[k] - x[k]);
            const double next = objective(trial, &trialDirection);
            if (std::isfinite(next) && next <= value + 1e-4 * std::min(slope, 0.0)) {
                value = next;
                accepted = true;
                break;
            }
        }
        if (!accepted) {
            if (sHistory.empty()) {
                stopReason = "line search failed after " + std::to_string(lineSearchAttempts) + " trials";
                break;
            }
            ++lineSearchRetries;
            sHistory.clear();   // retry once from the gradient
            yHistory.clear();
            continue;
        }

        const std::vector<double> gTrial = pack(trialDirection, true);
        std::vector<double> s(dim), y(dim);
        for (size_t k = 0; k < dim; ++k) { s[k] = xTrial[k] - x[k]; y[k] = gTrial[k] - g[k]; }
        const double sy = dot(s, y);
        if (sy > 1e-10 * std::sqrt(dot(s, s) * dot(y, y))) {
            sHistory.push_back(std::move(s));
            yHistory.push_back(std::move(y));
            if (sHistory.size() > memory) { sHistory.pop_front(); yHistory.pop_front(); }
        }
        x = std::move(xTrial);
        g = gTrial;
        current = trial;
        direction = trialDirection;
        ++result.leafIterations;
        if (value < bestValue) { best = current; bestValue = value; }
        bestHistory.push_back(bestValue);
        const int acceptedSteps = result.leafIterations - initialAcceptedSteps;
        if (refine.progressCallback && shouldLogOptimizationProgress(acceptedSteps, 10)) {
            DaoProgressInfo info;
            info.iteration = iteration + 1;
            info.acceptedSteps = acceptedSteps;
            info.preStepGradInf = maxGradient;
            info.stepSize = step;
            info.lsIters = lineSearchAttempts;
            report("", current, info, false, false);
        }
        if (shouldLogOptimizationProgress(acceptedSteps, 50)) {
            if (acceptedSteps == 1) logDaoProgressHeader(refine.iterations, memory);
            double totalMu = 0.0;
            for (const auto& cp : current.controlPoints) totalMu += cp.mu;
            DaoProgressInfo info;
            info.iteration = iteration + 1;
            info.acceptedSteps = acceptedSteps;
            info.objective = value;
            info.preStepGradInf = maxGradient;
            info.stepSize = step;
            info.lsIters = lineSearchAttempts;
            info.totalMu = totalMu;
            info.elapsedSec = std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count();
            logDaoProgress(info);
        }
        if (bestHistory.size() > 40 && bestHistory[bestHistory.size() - 1 - 40] - bestValue < 1e-3 * bestValue) {
            stopReason = "best cost improved less than 0.1% over 40 accepted steps";
            break;
        }
    }

    // Deliverable positions: snap to the MLC resolution (still inside the speed and MU bounds).
    result.sequence = best;
    if (!project(result.sequence, false)) {
        Logger::warn("DAO feasibility repair failed; reverting to the stage seed");
        result.sequence = seed;
    }
    if (refine.finalStage) {
        const ArcSequenceResult repaired = result.sequence;
        if (!project(result.sequence, true) || !std::isfinite(evaluate(result.sequence, nullptr))) {
            Logger::warn("DAO leaf-position snapping failed; keeping the unsnapped repaired sequence");
            result.sequence = repaired;
        }
    }
    value = evaluate(result.sequence, nullptr);
    result.finalLeafObjective = value;
    ArcSequencer::computeTiming(result.sequence, limits);
    if (refine.finalStage || refine.progressCallback) {
        DaoProgressInfo info;
        info.iteration = attemptedIterations;
        info.acceptedSteps = result.leafIterations - initialAcceptedSteps;
        report("DAO final (after feasibility repair)", result.sequence, info, refine.finalStage, true);
    }
    std::ostringstream summary;
    summary << "DAO stopped: " << stopReason << "\n  Iter " << attemptedIterations << "/" << refine.iterations
            << " | Acc " << result.leafIterations - initialAcceptedSteps
            << " | Restarts " << directionRestarts << "/" << lineSearchRetries << " (dir/LS)"
            << " | Time " << std::setprecision(4)
            << std::chrono::duration<double>(std::chrono::steady_clock::now() - startTime).count()
            << " s\n  DoseCost " << std::scientific << result.initialLeafObjective << " -> " << value
            << " | MU " << std::defaultfloat << std::setprecision(5) << seedMu << " -> " << result.sequence.totalMU
            << "\n  Delivery " << result.sequence.totalTimeSec << " s | LeafMax "
            << result.sequence.maxLeafSpeed << " mm/s | Slowed " << result.sequence.slowedControlPoints << " CPs";
    Logger::info(summary.str());
}

std::shared_ptr<FluenceDoseModel> buildStageModel(const ArcPipelineResult& stage, const Plan& plan,
                                                  const PatientData& patient, const Grid& grid,
                                                  const ArcOptimizationConfig& config) {
    PencilBeamEngine engine;
    engine.setOptions(config.engineOptions);
    FluenceModelOptions options;
    options.pixelMm = config.fluencePixelMm;
    options.backend = config.computeBackend;
    options.gpuMemoryBudgetMiB = config.gpuMemoryBudgetMiB;
    return std::make_shared<FluenceDoseModel>(
        engine.buildFluenceModel(plan, *stage.controlPointStf, patient, grid, options));
}

// Dose is linear in MU, so one dose calculation gives the global MU scale of the seed that minimizes the objectives.
void calibrateSeedMu(ArcPipelineResult& stage, const FluenceDoseOperator& doseOperator,
                     const BuiltObjectives& objectives, const ArcOptimizationConfig& config) {
    Logger::info("DAO seed: calibrating uniform MU against the dose objectives");
    for (auto& cp : stage.sequence.controlPoints) cp.mu = cp.aperture.weight = 1.0;
    const auto unitDose = doseOperator.computeDose(
        ArcPipeline::bixelWeights(*stage.controlPointStf, stage.sequence, config.integrateLeafMotion));
    auto cost = [&](double scale) {
        std::vector<double> dose(unitDose);
        for (double& value : dose) value *= scale;
        double total = 0.0;
        for (const auto* objective : objectives.ptrs) total += objective->compute(dose);
        return total;
    };
    // One-sided objectives (minimum dose) are flat above their optimum: ties go to the smaller scale.
    const double golden = 0.5 * (std::sqrt(5.0) - 1.0);
    double low = std::log(1e-2), high = std::log(1e7);
    double a = high - golden * (high - low), b = low + golden * (high - low);
    double fa = cost(std::exp(a)), fb = cost(std::exp(b));
    for (int i = 0; i < 80; ++i) {
        if (fa <= fb) { high = b; b = a; fb = fa; a = high - golden * (high - low); fa = cost(std::exp(a)); }
        else { low = a; a = b; fa = fb; b = low + golden * (high - low); fb = cost(std::exp(b)); }
    }
    const double scale = std::exp(0.5 * (low + high));
    for (auto& cp : stage.sequence.controlPoints) cp.mu = cp.aperture.weight = scale;
    Logger::info("DAO seed: uniform MU per control point " + std::to_string(scale));
}

// Progressive resolution direct aperture optimization (as in Eclipse's VMAT optimizer). The seed opens every leaf
// pair over the target projection of the sector beams (`sectorStf`) and its MU is calibrated against the objectives.
// The normal tissue objective is added once the first dose model is known (its dose rows decide which body voxels
// it can see).
void runProgressive(ArcPipelineResult& result, const Stf& sectorStf,
                    const Plan& plan, const PatientData& patient, const Grid& grid,
                    BuiltObjectives& objectives, const ArcOptimizationConfig& config,
                    const std::function<void(const std::string&)>& say,
                    const ArcProgressCallback& progressCallback) {
    const VmatSettings finalSettings = plan.getVmatSettings();
    const MachineConstraints limits = ArcPipeline::deliveryConstraints(plan);
    std::vector<double> spacings{finalSettings.controlPointSpacingDeg};
    while (spacings.back() * 2.0 <= config.coarsestSpacingDeg) spacings.push_back(spacings.back() * 2.0);
    std::reverse(spacings.begin(), spacings.end());

    ArcSequenceResult previous;
    double muLambda = 0.0;   // stage 1 runs without MU penalty; the penalty is derived at the first stage after it
    double smoothingReference = 0.0;   // likewise for the MU smoothing (0 = off)
    for (size_t level = 0; level < spacings.size(); ++level) {
        VmatSettings settings = finalSettings;
        settings.controlPointSpacingDeg = spacings[level];
        const bool first = level == 0, last = level + 1 == spacings.size();

        ArcPipelineResult stage;
        say("DAO stage " + std::to_string(level + 1) + "/" + std::to_string(spacings.size()) +
            (first ? ": creating conformal seed" : ": interpolating previous apertures and MU"));
        if (first) {
            stage.sequence = ArcSequencer::conformalSeed(sectorStf, settings, plan.getMachine(), config.seed);
        } else {
            std::vector<ArcControlPoint> geometry = ArcGeometry::controlPoints(settings);
            stage.sequence = ArcSequencer::resample(previous, geometry, limits);
        }
        say("resolution " + std::to_string(level + 1) + "/" + std::to_string(spacings.size()) + ": " +
            std::to_string(stage.sequence.controlPoints.size()) + " control points (" +
            std::to_string(spacings[level]) + " deg)");

        std::vector<ArcControlPoint> geometry;
        for (const auto& cp : stage.sequence.controlPoints) geometry.push_back(cp.geometry);
        say("building stage control point beams and preparing dose model");
        stage.controlPointStf = ArcPipeline::buildControlPointStf(plan, patient, geometry);

        std::shared_ptr<FluenceDoseModel> preparedModel;
        if (first) {
            preparedModel = buildStageModel(stage, plan, patient, grid, config);
            FluenceDoseOperator doseOperator(preparedModel, stage.controlPointStf);
            ObjectiveBuilder::addNormalTissueObjective(objectives, config.base.nto, config.base.targetDose,
                                                       patient, grid, &doseOperator);
            if (objectives.ptrs.empty()) throw std::invalid_argument("ArcPipeline: DAO needs at least one objective");
            logObjectiveConfiguration("ArcPipeline", objectives.ptrs);
            calibrateSeedMu(stage, doseOperator, objectives, config);
        }

        RefineOptions refine;
        refine.stage = static_cast<int>(level + 1);
        refine.stages = static_cast<int>(spacings.size());
        refine.progressCallback = progressCallback;
        refine.iterations = last ? 2 * config.iterationsPerResolution : config.iterationsPerResolution;
        refine.reportInitial = first;
        refine.finalStage = last;
        refine.muPenalty = config.muPenalty;
        const double width = stage.sequence.controlPoints.front().geometry.width;
        const double gantryMax = limits.gantryRotationSpeed[1] > 0.0 ? limits.gantryRotationSpeed[1] : 6.0;
        refine.maxMoveMm = std::max(4.0 * config.leafStepMm, 0.25 * limits.leafSpeed[1] * width / gantryMax);
        // Stage 1 finds the MU scale freely; the penalty is then fixed from the first plan-like stage.
        const bool derivesLambda = level == (spacings.size() > 1 ? 1 : 0);
        double lambda = derivesLambda ? -1.0 : muLambda;
        refine.muLambda = &lambda;
        if (derivesLambda) smoothingReference = -1.0;
        refine.smoothingReference = &smoothingReference;
        refineLeaves(stage, plan, patient, grid, objectives, config, preparedModel, refine);
        result.progressHistory.insert(result.progressHistory.end(),
            stage.progressHistory.begin(), stage.progressHistory.end());
        if (derivesLambda) muLambda = lambda;

        if (first) result.initialLeafObjective = stage.initialLeafObjective;
        result.finalLeafObjective = stage.finalLeafObjective;
        result.leafIterations += stage.leafIterations;
        previous = stage.sequence;
        if (last) {
            result.sequence = std::move(stage.sequence);
            result.controlPointStf = std::move(stage.controlPointStf);
        }
    }
}

} // namespace

std::shared_ptr<Stf> ArcPipeline::buildControlPointStf(const Plan& plan,
                                                       const PatientData& patientData,
                                                       const std::vector<ArcControlPoint>& controlPoints) {
    const auto& props = plan.getStfProperties();
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    if (!props.isoCenters.empty()) iso = props.isoCenters[0];

    std::vector<double> gantry, couch;
    for (const auto& cp : controlPoints) {
        gantry.push_back(cp.gantryAngle);
        couch.push_back(cp.couchAngle);
    }

    PhotonIMRTStfGenerator generator(0.0, 360.0, 360.0, props.bixelWidth, iso);
    generator.setMachine(plan.getMachine());
    generator.setRadiationMode(plan.getRadiationMode());
    generator.setGantryAngles(gantry);
    generator.setCouchAngles(couch);

    if (patientData.hasValidCT() && patientData.hasStructures()) {
        const auto& grid = patientData.getCTVolume()->getGrid();
        generator.setGrid(grid);
        generator.setStructureSet(*patientData.getStructureSet());
        generator.setCTResolution(grid.getSpacing());
    }
    return std::make_shared<Stf>(generator.generateStf());
}

std::vector<double> ArcPipeline::bixelWeights(const Stf& cpStf, const ArcSequenceResult& sequence,
                                              bool integrateLeafMotion) {
    if (cpStf.getCount() != sequence.controlPoints.size())
        throw std::invalid_argument("ArcPipeline::bixelWeights: STF and control point counts differ");

    std::vector<double> weights(cpStf.getTotalNumOfBixels(), 0.0);
    size_t offset = 0;
    for (size_t j = 0; j < cpStf.getCount(); ++j) {
        const Beam& beam = *cpStf.getBeam(j);
        const auto& cp = sequence.controlPoints[j];
        if (cp.mu > 0.0) {
            const auto open = integrateLeafMotion ?
                ApertureFluence::integratedOpenFractions(beam, neighbour(sequence, j, -1), cp.aperture,
                    neighbour(sequence, j, 1), cp.leafPairBoundariesZ) :
                ApertureFluence::bixelOpenFractions(beam, cp.aperture, cp.leafPairBoundariesZ);
            for (size_t r = 0; r < open.size(); ++r) weights[offset + r] = open[r] * cp.mu;
        }
        offset += beam.getTotalNumOfBixels();
    }
    return weights;
}

double ArcPipeline::machineMuPerWeight(const Plan& plan) { return plan.machineMuPerWeight(); }

MachineConstraints ArcPipeline::deliveryConstraints(const Plan& plan) {
    MachineConstraints limits = plan.getMachine().getConstraints();
    const double muPerWeight = machineMuPerWeight(plan);
    for (double& rate : limits.monitorUnitRate) rate /= muPerWeight;
    return limits;
}

ArcPipelineResult ArcPipeline::run(const Stf& sectorStf,
                                   const Plan& plan,
                                   const PatientData& patientData,
                                   const Grid& doseGrid,
                                   BuiltObjectives objectives,
                                   const ArcOptimizationConfig& config,
                                   std::function<void(const std::string&)> stage,
                                   ArcProgressCallback progressCallback) {
    auto say = [&](const std::string& s) { Logger::info("ArcPipeline: " + s); if (stage) stage(s); };
    if (!plan.isVmat()) throw std::invalid_argument("ArcPipeline: the plan is not a VMAT plan");
    if (config.base.targetDose <= 0.0 || !std::isfinite(config.base.targetDose))
        throw std::invalid_argument("ArcPipeline: invalid prescription");

    ArcPipelineResult result;
    say("starting progressive DAO, dose voxels=" + std::to_string(doseGrid.getNumVoxels()) +
        ", leaf-motion integration=" + (config.integrateLeafMotion ? "enabled" : "disabled"));
    auto engineOptions = config.engineOptions;
    engineOptions.gpuRayTracing = config.computeBackend != DoseComputeBackend::CPU;
    engineOptions.requireGpuRayTracing = config.computeBackend == DoseComputeBackend::CUDA;
#ifndef OPTIRAD_HAS_CUDA
    if (engineOptions.requireGpuRayTracing)
        throw std::runtime_error("CUDA arc dose requested but CUDA support was not built");
#endif

    runProgressive(result, sectorStf, plan, patientData, doseGrid, objectives, config, say, progressCallback);
    auto engine = DoseEngineFactory::create("PencilBeam");
    engine->setOptions(engineOptions);
    say("computing delivered dose");
    const std::vector<double> weights = bixelWeights(*result.controlPointStf, result.sequence,
                                                     config.integrateLeafMotion);
    result.dose = std::make_shared<DoseMatrix>(
        engine->calculateDoseDirect(plan, *result.controlPointStf, patientData, doseGrid, weights));
    result.stats = PlanAnalysis::computeStats(*result.dose, patientData, doseGrid, config.base.targetDose);
    say("complete: control points=" + std::to_string(result.sequence.controlPoints.size()) +
        ", MU=" + std::to_string(result.sequence.totalMU) + " (~" +
        std::to_string(result.sequence.totalMU * machineMuPerWeight(plan)) + " machine MU/fraction)" +
        ", delivery time=" + std::to_string(result.sequence.totalTimeSec) + " s");
    return result;
}

} // namespace optirad
