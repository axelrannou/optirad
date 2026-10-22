#include "ArcPipeline.hpp"

#include "dose/DoseEngineFactory.hpp"
#include "dose/FluenceDoseOperator.hpp"
#include "dose/engines/PencilBeamEngine.hpp"
#include "optimization/OptimizerFactory.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "sequencing/ApertureFluence.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "utils/Logger.hpp"

#include <stdexcept>
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>

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
    int iterations = 0;
    double muCapFactor = 0.0;       // MU of a control point may grow to this multiple of its seed; 0 = rate limit only
    double maxMoveMm = 0.0;         // largest leaf tip move per iteration; 0 = 4 x leafStepMm
    double* muLambda = nullptr;     // objective units per MU: < 0 on entry = derive from the initial objective
    double muPenalty = 0.0;         // fraction of the starting objective that the total MU costs
    double hotspotThresholdPct = 1.01;   // fraction of the prescription
    double hotspotPenalty = 50000.0;
    bool finalStage = true;         // snap to the MLC resolution and report timing at the end
};

void refineLeaves(ArcPipelineResult& result, const Plan& plan, const PatientData& patient,
                  const Grid& grid, const BuiltObjectives& objectives, const ArcOptimizationConfig& config,
                  std::shared_ptr<FluenceDoseModel> preparedModel, const RefineOptions& refine) {
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
    const auto& limits = plan.getMachine().getConstraints();
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
            const double maxMU = std::min(limits.monitorUnitRate[1] * duration,
                refine.muCapFactor > 0.0 ? std::max(original.mu, 1e-6) * refine.muCapFactor :
                    std::numeric_limits<double>::infinity());
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
    auto evaluate = [&](const ArcSequenceResult& sequence, std::vector<double>* gradient) {
        const auto weights = ArcPipeline::bixelWeights(*result.controlPointStf, sequence, config.integrateLeafMotion);
        double value = 0.0;
        std::vector<double> unusedGradient;
        if (!supported || !doseOperator.evaluateObjectives(weights, descriptions, config.base.targetDose,
            refine.hotspotThresholdPct, config.base.ntoEnabled ? refine.hotspotPenalty : 0.0,
            value, gradient ? *gradient : unusedGradient, true, gradient != nullptr)) {
        const auto dose = doseOperator.computeDose(weights);
        std::vector<double> voxelGradient(gradient ? dose.size() : 0, 0.0);
        for (const auto* objective : objectives.ptrs) {
            value += objective->compute(dose);
            if (gradient) {
                const auto contribution = objective->gradient(dose);
                for (size_t voxel = 0; voxel < dose.size(); ++voxel) voxelGradient[voxel] += contribution[voxel];
            }
        }
        if (config.base.ntoEnabled) {
            const double prescription = config.base.targetDose;
            for (size_t voxel = 0; voxel < dose.size(); ++voxel) {
                const double excess = dose[voxel] - prescription * refine.hotspotThresholdPct;
                if (excess <= 0.0) continue;
                const double normalized = excess / prescription;
                value += refine.hotspotPenalty * normalized * normalized * (1.0 + 10.0 * normalized) * excess * excess;
                if (gradient) voxelGradient[voxel] += refine.hotspotPenalty * excess * normalized * normalized *
                    (4.0 + 50.0 * normalized);
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
            value += config.muSmoothingWeight * difference * difference;
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
        muHigh[index] = std::min(limits.monitorUnitRate[1] * original.timeSec,
            refine.muCapFactor > 0.0 ? std::max(original.mu, 1e-6) * refine.muCapFactor :
                std::numeric_limits<double>::infinity());
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
            const double derivative = 2.0 * config.muSmoothingWeight * (cp.mu - previous.mu);
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
            direction->controlPoints[index].mu = config.optimizeWeights ?
                direction->controlPoints[index].mu * muScale[index] : 0.0;
        return value;
    };

    result.initialLeafObjective = evaluate(result.sequence, nullptr);
    if (deriveLambda && config.optimizeWeights && refine.muPenalty > 0.0 && seedMu > 0.0)
        muLambda = refine.muPenalty * result.initialLeafObjective / seedMu;

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
    for (int iteration = 0; iteration < refine.iterations; ++iteration) {
        double maxGradient = 0.0;
        for (double v : g) maxGradient = std::max(maxGradient, std::abs(v));
        if (maxGradient <= 0.0) break;

        std::vector<double> d(g);
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
        freeze(d);
        if (dot(d, g) >= 0.0) {   // not a descent direction: restart from the gradient
            sHistory.clear();
            yHistory.clear();
            for (size_t k = 0; k < dim; ++k) d[k] = -g[k] * config.leafStepMm / maxGradient;
            freeze(d);
            if (dot(d, g) >= 0.0) break;   // stationary on the feasible box
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
        std::vector<double> xTrial(dim);
        for (int attempt = 0; attempt < 20; ++attempt, step *= 0.5) {
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
            if (sHistory.empty()) break;
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
        if (iteration >= 30 && bestHistory[bestHistory.size() - 1 - 15] - bestValue < 1e-3 * bestValue) break;
    }

    // Deliverable positions: snap to the MLC resolution (still inside the speed and MU bounds).
    result.sequence = best;
    if (!project(result.sequence, false)) result.sequence = seed;   // the speed penalty is soft: repair exactly
    if (refine.finalStage) {
        const ArcSequenceResult repaired = result.sequence;
        if (!project(result.sequence, true) || !std::isfinite(evaluate(result.sequence, nullptr)))
            result.sequence = repaired;
    }
    value = evaluate(result.sequence, nullptr);
    result.finalLeafObjective = value;
    ArcSequencer::computeTiming(result.sequence, limits);
    Logger::info("Leaf refinement: " + std::to_string(result.leafIterations) + " accepted steps, objective " +
        std::to_string(result.initialLeafObjective) + " -> " + std::to_string(value) + ", MU " +
        std::to_string(seedMu) + " -> " + std::to_string(result.sequence.totalMU));
}

// Progressive resolution direct aperture optimization (see ArcStrategy::ProgressiveDao).
void runProgressive(ArcPipelineResult& result, const std::vector<double>& fmoWeights, const Stf& fmoStf,
                    const Plan& plan, const PatientData& patient, const Grid& grid,
                    const BuiltObjectives& objectives, const ArcOptimizationConfig& config,
                    const std::function<void(const std::string&)>& say) {
    const VmatSettings finalSettings = plan.getVmatSettings();
    const auto& limits = plan.getMachine().getConstraints();
    std::vector<double> spacings{finalSettings.controlPointSpacingDeg};
    while (spacings.back() * 2.0 <= config.daoCoarsestSpacingDeg) spacings.push_back(spacings.back() * 2.0);
    std::reverse(spacings.begin(), spacings.end());

    ArcSequenceResult previous;
    double muLambda = 0.0;   // stage 1 runs without MU penalty; the penalty is derived at the first stage after it
    for (size_t level = 0; level < spacings.size(); ++level) {
        VmatSettings settings = finalSettings;
        settings.controlPointSpacingDeg = spacings[level];
        const bool first = level == 0, last = level + 1 == spacings.size();

        ArcPipelineResult stage;
        if (first) {
            ArcSequencerOptions seed = config.sequencer;
            seed.method = ArcSequencerOptions::Method::Conformal;
            stage.sequence = ArcSequencer::sequence(fmoWeights, fmoStf, settings, plan.getMachine(), seed);
        } else {
            std::vector<ArcControlPoint> geometry = ArcGeometry::controlPoints(settings);
            stage.sequence = ArcSequencer::resample(previous, geometry, limits);
        }
        say("resolution " + std::to_string(level + 1) + "/" + std::to_string(spacings.size()) + ": " +
            std::to_string(stage.sequence.controlPoints.size()) + " control points (" +
            std::to_string(spacings[level]) + " deg)");

        std::vector<ArcControlPoint> geometry;
        for (const auto& cp : stage.sequence.controlPoints) geometry.push_back(cp.geometry);
        stage.controlPointStf = ArcPipeline::buildControlPointStf(plan, patient, geometry);

        RefineOptions refine;
        refine.iterations = last ? 2 * config.daoIterations : config.daoIterations;
        refine.finalStage = last;
        refine.muPenalty = config.daoMuPenalty;
        // Continuation: the hot-spot term starts like the fluence optimization's and tightens to the arc setting.
        const double t = spacings.size() > 1 ? static_cast<double>(level) / static_cast<double>(spacings.size() - 1) : 1.0;
        refine.hotspotThresholdPct = (1.0 - t) * config.base.ntoThresholdPct + t * config.daoHotspotThresholdPct;
        refine.hotspotPenalty = std::exp((1.0 - t) * std::log(std::max(config.base.ntoPenalty, 1.0)) +
                                         t * std::log(std::max(config.daoHotspotPenalty, 1.0)));
        const double width = stage.sequence.controlPoints.front().geometry.width;
        const double gantryMax = limits.gantryRotationSpeed[1] > 0.0 ? limits.gantryRotationSpeed[1] : 6.0;
        refine.maxMoveMm = std::max(4.0 * config.leafStepMm, 0.25 * limits.leafSpeed[1] * width / gantryMax);
        // Stage 1 finds the MU scale freely; the penalty is then fixed from the first plan-like stage.
        const bool derivesLambda = level == (spacings.size() > 1 ? 1 : 0);
        double lambda = derivesLambda ? -1.0 : muLambda;
        refine.muLambda = &lambda;
        refineLeaves(stage, plan, patient, grid, objectives, config, nullptr, refine);
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

ArcPipelineResult ArcPipeline::run(const std::vector<double>& fmoWeights,
                                   const Stf& fmoStf,
                                   const Plan& plan,
                                   const PatientData& patientData,
                                   const Grid& doseGrid,
                                   BuiltObjectives objectives,
                                   const ArcOptimizationConfig& config,
                                   IterationCallback iterCallback,
                                   std::function<void(const std::string&)> stage) {
    auto say = [&](const std::string& s) { Logger::info("ArcPipeline: " + s); if (stage) stage(s); };
    if (!plan.isVmat()) throw std::invalid_argument("ArcPipeline: the plan is not a VMAT plan");
    if (config.maxWeightFactor < 1.0 || !std::isfinite(config.maxWeightFactor) ||
        config.base.targetDose <= 0.0 || !std::isfinite(config.base.targetDose))
        throw std::invalid_argument("ArcPipeline: invalid MU growth factor or prescription");

    ArcPipelineResult result;

    if (config.strategy == ArcStrategy::ProgressiveDao) {
        if (objectives.ptrs.empty() || !config.optimizeLeaves)
            throw std::invalid_argument("ArcPipeline: progressive optimization needs objectives and leaf optimization");
        runProgressive(result, fmoWeights, fmoStf, plan, patientData, doseGrid, objectives, config, say);
        auto engine = DoseEngineFactory::create("PencilBeam");
        engine->setOptions(config.engineOptions);
        say("computing delivered dose");
        const std::vector<double> weights = bixelWeights(*result.controlPointStf, result.sequence,
                                                         config.integrateLeafMotion);
        result.dose = std::make_shared<DoseMatrix>(
            engine->calculateDoseDirect(plan, *result.controlPointStf, patientData, doseGrid, weights));
        result.stats = PlanAnalysis::computeStats(*result.dose, patientData, doseGrid, config.base.targetDose);
        return result;
    }

    say("arc sequencing");
    result.sequence = ArcSequencer::sequence(fmoWeights, fmoStf, plan.getVmatSettings(), plan.getMachine(),
                                             config.sequencer);

    say("building " + std::to_string(result.sequence.controlPoints.size()) + " control point beams");
    std::vector<ArcControlPoint> geometry;
    for (const auto& cp : result.sequence.controlPoints) geometry.push_back(cp.geometry);
    result.controlPointStf = buildControlPointStf(plan, patientData, geometry);

    auto engine = DoseEngineFactory::create("PencilBeam");
    engine->setOptions(config.engineOptions);
    std::shared_ptr<FluenceDoseModel> apertureModel;

    if (config.optimizeWeights && !objectives.ptrs.empty()) {
        const std::vector<double> baseWeights = bixelWeights(*result.controlPointStf, result.sequence,
            config.integrateLeafMotion);
        std::shared_ptr<IDoseOperator> apertureDose;
#ifdef OPTIRAD_HAS_CUDA
        if (config.computeBackend != DoseComputeBackend::CPU && config.engineOptions.divergentBixels) {
            say("preparing CUDA matrix-free aperture dose model");
            PencilBeamEngine modelEngine;
            modelEngine.setOptions(config.engineOptions);
            FluenceModelOptions options;
            options.pixelMm = config.fluencePixelMm;
            options.backend = config.computeBackend;
            options.gpuMemoryBudgetMiB = config.gpuMemoryBudgetMiB;
            apertureModel = std::make_shared<FluenceDoseModel>(modelEngine.buildFluenceModel(
                plan, *result.controlPointStf, patientData, doseGrid, options));
            if (apertureModel->usesCuda()) {
                auto bixelOperator = std::make_shared<FluenceDoseOperator>(apertureModel, result.controlPointStf);
                apertureDose = std::make_shared<ApertureDoseOperator>(bixelOperator, *result.controlPointStf, baseWeights);
                auto doseOptions = config.engineOptions;
                doseOptions.gpuRayTracing = true;
                engine->setOptions(doseOptions);
            } else apertureModel.reset();
        }
#else
        if (config.computeBackend == DoseComputeBackend::CUDA)
            throw std::runtime_error("CUDA aperture optimization requested but CUDA support was not built");
#endif
        if (!apertureDose) {
            say("computing aperture dose columns (CPU operator)");
            apertureDose = std::make_shared<DoseInfluenceMatrix>(
                engine->calculateBeamDij(plan, *result.controlPointStf, patientData, doseGrid, baseWeights));
        }

        const size_t n = result.sequence.controlPoints.size();
        std::vector<double> x0(n, 0.0);
        for (size_t j = 0; j < n; ++j)
            if (result.sequence.controlPoints[j].mu > 0.0) x0[j] = 1.0;

        // Neighbours: adjacent control points of the same arc.
        std::vector<std::vector<int>> neighbors(n);
        for (size_t j = 0; j + 1 < n; ++j) {
            const auto& a = result.sequence.controlPoints[j].geometry;
            const auto& b = result.sequence.controlPoints[j + 1].geometry;
            if (a.arcIndex == b.arcIndex) {
                neighbors[j].push_back(static_cast<int>(j + 1));
                neighbors[j + 1].push_back(static_cast<int>(j));
            }
        }

        say(apertureModel ? "optimizing aperture weights (CUDA dose/objectives/gradient)" :
            "optimizing aperture weights (CPU stored columns)");
        auto optimizer = OptimizerFactory::create("LBFGS");
        optimizer->setMaxIterations(config.weightIterations);
        optimizer->setTolerance(config.base.tolerance);
        if (iterCallback) optimizer->setIterationCallback(std::move(iterCallback));
        if (auto* lbfgs = dynamic_cast<LBFGSOptimizer*>(optimizer.get())) {
            lbfgs->setMaxFluence(config.maxWeightFactor);
            lbfgs->setInitialWeights(x0);
            lbfgs->setPrescriptionDose(config.base.targetDose);
            if (config.base.ntoEnabled) {
                lbfgs->setHotspotThreshold(config.hotspotThresholdPct);
                lbfgs->setHotspotPenalty(config.hotspotPenalty);
                lbfgs->setConsistentHotspotGradient(true);
            }
            if (config.muSmoothingWeight > 0.0) lbfgs->setSpatialSmoothing(config.muSmoothingWeight, neighbors);
            if (config.muPenalty > 0.0) lbfgs->setRelativeL1Regularization(config.muPenalty);
        }

        std::vector<Constraint> constraints;
        OptimizationResult opt = optimizer->optimize(*apertureDose, objectives.ptrs, constraints);
        result.iterations = opt.iterations;
        result.finalObjective = opt.finalObjective;
        result.converged = opt.converged;

        for (size_t j = 0; j < n; ++j) {
            auto& cp = result.sequence.controlPoints[j];
            cp.mu *= opt.weights[j];
            cp.aperture.weight = cp.mu;
        }
        ArcSequencer::computeTiming(result.sequence, plan.getMachine().getConstraints());
    }

    if (config.optimizeLeaves && !objectives.ptrs.empty()) {
        say("refining leaf positions and MU (matrix-free)");
        RefineOptions refine;
        refine.hotspotThresholdPct = config.hotspotThresholdPct;
        refine.hotspotPenalty = config.hotspotPenalty;
        refine.muPenalty = config.muPenalty;
        refine.iterations = config.leafIterations;
        refine.muCapFactor = config.maxWeightFactor;
        double lambda = -1.0;
        refine.muLambda = &lambda;
        refineLeaves(result, plan, patientData, doseGrid, objectives, config, apertureModel, refine);
    }

    say("computing delivered dose");
    const std::vector<double> finalWeights = bixelWeights(*result.controlPointStf, result.sequence,
        config.integrateLeafMotion);
    result.dose = std::make_shared<DoseMatrix>(
        engine->calculateDoseDirect(plan, *result.controlPointStf, patientData, doseGrid, finalWeights));
    result.stats = PlanAnalysis::computeStats(*result.dose, patientData, doseGrid, config.base.targetDose);
    return result;
}

} // namespace optirad
