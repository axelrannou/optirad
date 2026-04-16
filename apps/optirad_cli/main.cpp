#include "io/DicomImporter.hpp"
#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Stf.hpp"
#include "core/workflow/WorkflowState.hpp"
#include "core/workflow/PlanBuilder.hpp"
#include "core/workflow/DoseCalculationPipeline.hpp"
#include "core/workflow/PhaseSpaceBuilder.hpp"
#include "core/workflow/OptimizationPipeline.hpp"
#include "core/workflow/LeafSequencingPipeline.hpp"
#include "optimization/ObjectiveProtocol.hpp"
#include "dose/PlanAnalysis.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <cstdlib>
#include <chrono>
#include <algorithm>
#include <numeric>
#include <iomanip>

using namespace optirad;

// Forward declarations
int loadPhaseSpace(const std::vector<std::string>& args, WorkflowState& state);
int doseCalc(const std::vector<std::string>& args, WorkflowState& state);
int optimize(const std::vector<std::string>& args, WorkflowState& state);
int leafSeq(const std::vector<std::string>& args, WorkflowState& state);

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " <command> [options]\n\n"
              << "Commands:\n"
              << "  load <dicom_dir>                Load and inspect DICOM directory\n"
              << "  plan [options]                   Generate a treatment plan and STF (requires DICOM data first)\n"
              << "  doseCalc [options]               Calculate Dij (requires plan)\n"
              << "  optimize [options]               Run optimization (requires Dij)\n"
              << "  leafSeq [options]                Run leaf sequencing (requires optimization)\n"
              << "  loadPhaseSpace [options]          Load phase-space beam source (phase-space machines only)\n"
              << "  interactive                      Enter interactive mode\n"
              << "  help                             Show this help message\n\n"
              << "Plan options:\n"
              << "  --mode <photons|protons>         Radiation mode (default: photons)\n"
              << "  --machine <name>                 Machine name (default: Generic, Varian_TrueBeam6MV)\n"
              << "  --fractions <n>                  Number of fractions (default: 30)\n"
              << "  --gantry-start <deg>             Gantry start angle (default: 0)\n"
              << "  --gantry-step <deg>              Gantry angle step (default: 4)\n"
              << "  --gantry-stop <deg>              Gantry stop angle exclusive (default: 360)\n"
              << "  --gantry-angles <a1 a2 ...>      Explicit gantry angle list (space-separated)\n"
              << "  --couch-start <deg>              Couch start angle (default: 0)\n"
              << "  --couch-step <deg>               Couch angle step (default: 0 = single angle)\n"
              << "  --couch-stop <deg>               Couch stop angle exclusive (default: 0)\n"
              << "  --couch-angles <a1 a2 ...>       Explicit couch angle list (space-separated)\n"
              << "  --bixel-width <mm>               Bixel width (default: 7)\n\n"
              << "Dose calc options:\n"
              << "  --dose-resolution <mm>           Dose grid resolution (default: 2.5)\n"
              << "  --no-cache                       Disable Dij cache\n\n"
              << "Optimize options:\n"
              << "  --max-iter <n>                   Max iterations (default: 500)\n"
              << "  --tolerance <val>                Convergence tolerance (default: 1e-5)\n"
              << "  --target-dose <Gy>               Prescribed dose for targets (default: 60)\n"
              << "  --oar-max-dose <Gy>              Max dose for OARs (default: 30)\n\n"
              << "Leaf sequencing options:\n"
              << "  --num-levels <n>                 Intensity quantization levels (default: 15)\n"
              << "  --min-segment-mu <MU>            Min MU per segment (default: 0)\n\n"
              << "Phase-space options:\n"
              << "  --collimator <deg>               Collimator angle (default: 0)\n"
              << "  --max-particles <n>              Max particles per beam (default: 1000000)\n"
              << "  --viz-sample <n>                 Visualization sample per beam (default: 100000)\n"
              << "  (Couch angles come from plan's --couch-start/step/stop or --couch-angles)\n";
}

int loadDicom(const std::string& path, WorkflowState& state) {
    std::cout << "=== OptiRad DICOM Loader ===\n\n";
    std::cout << "Loading DICOM from: " << path << "\n\n";

    DicomImporter importer;

    if (!importer.canImport(path)) {
        std::cerr << "Error: Cannot import from path: " << path << "\n";
        return 1;
    }

    auto patientData = importer.importAll(path);

    if (!patientData) {
        std::cerr << "Error: Failed to import patient data\n";
        return 1;
    }

    // Store in shared state - move unique_ptr into shared_ptr
    state.patientData = std::move(patientData);

    // Invalidate derived state to avoid stale reuse across patients
    state.resetPlan();

    // Display patient info
    if (auto* patient = state.patientData->getPatient()) {
        std::cout << "\nPatient Information:\n";
        std::cout << "  Name: " << patient->getName() << "\n";
        std::cout << "  ID:   " << patient->getID() << "\n";
    }

    // Display CT volume info
    if (auto* ct = state.patientData->getCTVolume()) {
        const auto& grid = ct->getGrid();
        auto dims = grid.getDimensions();
        auto spacing = grid.getSpacing();
        auto origin = grid.getOrigin();
        auto orientation = grid.getImageOrientation();

        std::cout << "\nCT Volume:\n";
        std::cout << "  Dimensions: " << dims[0] << " x " << dims[1] << " x " << dims[2] << "\n";
        std::cout << "  Spacing:    " << spacing[0] << " x " << spacing[1] << " x " << spacing[2] << " mm\n";
        std::cout << "  Origin:     " << origin[0] << " x " << origin[1] << " x " << origin[2] << " mm\n";
        std::cout << "  Voxels:     " << ct->size() << "\n";

        // Get coordinate arrays
        auto x_coords = ct->getXCoordinates();
        auto y_coords = ct->getYCoordinates();
        auto z_coords = ct->getZCoordinates();
        
        // Print first 5 and last element
        std::cout << "              x: [";
        for (size_t i = 0; i < std::min(size_t(5), x_coords.size()); ++i) {
            std::cout << x_coords[i];
            if (i < 4 && i < x_coords.size() - 1) std::cout << " ";
        }
        if (x_coords.size() > 5) {
            std::cout << " ... " << x_coords.back();
        }
        std::cout << "] (1×" << x_coords.size() << " double)\n";
        
        std::cout << "              y: [";
        for (size_t i = 0; i < std::min(size_t(5), y_coords.size()); ++i) {
            std::cout << y_coords[i];
            if (i < 4 && i < y_coords.size() - 1) std::cout << " ";
        }
        if (y_coords.size() > 5) {
            std::cout << " ... " << y_coords.back();
        }
        std::cout << "] (1×" << y_coords.size() << " double)\n";
        
        std::cout << "              z: [";
        for (size_t i = 0; i < std::min(size_t(5), z_coords.size()); ++i) {
            std::cout << z_coords[i];
            if (i < 4 && i < z_coords.size() - 1) std::cout << " ";
        }
        if (z_coords.size() > 5) {
            std::cout << " ... " << z_coords.back();
        }
        std::cout << "] (1×" << z_coords.size() << " double)\n";

        int16_t minHU = 32767, maxHU = -32768;
        for (size_t i = 0; i < ct->size(); ++i) {
            minHU = std::min(minHU, ct->data()[i]);
            maxHU = std::max(maxHU, ct->data()[i]);
        }
        std::cout << "  HU Range:   [" << minHU << ", " << maxHU << "]\n";
    }

    // Display structures
    if (auto* structures = state.patientData->getStructureSet()) {
        std::cout << "\nStructures: " << structures->getCount() << "\n";
        for (size_t i = 0; i < structures->getCount(); ++i) {
            const auto* s = structures->getStructure(i);
            if (s) {
                auto color = s->getColor();
                std::cout << "  - " << s->getName()
                         << " (" << s->getType() << ")"
                         << " [" << s->getContourCount() << " contours]"
                         << " RGB(" << (int)color[0] << "," << (int)color[1] << "," << (int)color[2] << ")\n";
            }
        }
    }

    std::cout << "\n=== DICOM import complete ===\n";
    return 0;
}

int createPlan(const std::vector<std::string>& args, WorkflowState& state) {
    if (!state.patientData) {
        std::cerr << "Error: No patient data loaded. Use 'load <dicom_dir>' first.\n";
        return 1;
    }

    // Parse into PlanConfig
    PlanConfig config;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--mode" && i + 1 < args.size()) {
            config.radiationMode = args[++i];
        } else if (args[i] == "--machine" && i + 1 < args.size()) {
            config.machineName = args[++i];
        } else if (args[i] == "--fractions" && i + 1 < args.size()) {
            config.numFractions = std::stoi(args[++i]);
        } else if (args[i] == "--gantry-start" && i + 1 < args.size()) {
            config.gantryStart = std::stod(args[++i]);
        } else if (args[i] == "--gantry-step" && i + 1 < args.size()) {
            config.gantryStep = std::stod(args[++i]);
        } else if (args[i] == "--gantry-stop" && i + 1 < args.size()) {
            config.gantryStop = std::stod(args[++i]);
        } else if (args[i] == "--couch-start" && i + 1 < args.size()) {
            config.couchStart = std::stod(args[++i]);
        } else if (args[i] == "--couch-step" && i + 1 < args.size()) {
            config.couchStep = std::stod(args[++i]);
        } else if (args[i] == "--couch-stop" && i + 1 < args.size()) {
            config.couchStop = std::stod(args[++i]);
        } else if (args[i] == "--gantry-angles") {
            while (i + 1 < args.size() && args[i + 1][0] != '-') {
                config.gantryAngles.push_back(std::stod(args[++i]));
            }
        } else if (args[i] == "--couch-angles") {
            while (i + 1 < args.size() && args[i + 1][0] != '-') {
                config.couchAngles.push_back(std::stod(args[++i]));
            }
        } else if (args[i] == "--bixel-width" && i + 1 < args.size()) {
            config.bixelWidth = std::stod(args[++i]);
        }
    }

    try {
        auto result = PlanBuilder::build(config, state.patientData);
        state.plan = result.plan;
        state.stfProps = result.stfProps;
        state.stf = result.stf;
        state.resetDij();

        state.plan->printSummary();
        if (state.stf) {
            state.stf->printSummary();
        } else {
            std::cout << "\nPhase-space machine: use 'loadPhaseSpace' to build beam sources.\n";
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

int runInteractive(WorkflowState& state) {
    std::cout << "=== OptiRad Interactive Mode ===\n";
    std::cout << "Type 'help' for commands, 'quit' to exit.\n\n";

    std::string line;
    while (true) {
        std::cout << "optirad> ";
        if (!std::getline(std::cin, line)) break;

        // Tokenize
        std::istringstream iss(line);
        std::vector<std::string> tokens;
        std::string token;
        while (iss >> token) tokens.push_back(token);

        if (tokens.empty()) continue;

        const auto& cmd = tokens[0];

        if (cmd == "quit" || cmd == "exit") {
            break;
        } else if (cmd == "help") {
            std::cout << "Commands:\n"
                      << "  load <dicom_dir>              Load DICOM data\n"
                      << "  plan [options]                Create treatment plan and generate STF (requires DICOM data)\n"
                      << "  plan-help                     Show plan options and examples\n"
                      << "  doseCalc [options]            Calculate Dij (requires plan)\n"
                      << "  optimize [options]            Run optimization (requires Dij)\n"
                      << "  leafSeq [options]             Run leaf sequencing (requires optimization)\n"
                      << "  loadPhaseSpace [options]       Load phase-space beam source (PSF machines only)\n"
                      << "  phsp-info                     Show phase-space source metrics\n"
                      << "  info                          Show current state and available STF properties\n"
                      << "  quit                          Exit\n";
        } else if (cmd == "load") {
            if (tokens.size() < 2) {
                std::cerr << "Error: Missing DICOM directory path\n";
            } else {
                loadDicom(tokens[1], state);
            }
        } else if (cmd == "plan") {
            std::vector<std::string> planArgs(tokens.begin() + 1, tokens.end());
            createPlan(planArgs, state);
        } else if (cmd == "plan-help") {
            std::cout << "\n=== Plan Command Options ===\n"
                      << "Usage: plan [options]\n\n"
                      << "Options:\n"
                      << "  --mode <photons|protons>      Radiation mode (default: photons)\n"
                      << "  --machine <name>              Machine name (default: Generic)\n"
                      << "  --fractions <n>               Number of fractions (default: 30)\n"
                      << "  --gantry-start <deg>          Gantry start angle (default: 0)\n"
                      << "  --gantry-step <deg>           Gantry angle step (default: 4)\n"
                      << "  --gantry-stop <deg>           Gantry stop angle exclusive (default: 360)\n"
                      << "  --gantry-angles <a1 a2 ...>   Explicit gantry angle list (space-separated)\n"
                      << "  --couch-start <deg>           Couch start angle (default: 0)\n"
                      << "  --couch-step <deg>            Couch angle step (default: 0 = single angle)\n"
                      << "  --couch-stop <deg>            Couch stop angle exclusive (default: 0)\n"
                      << "  --couch-angles <a1 a2 ...>    Explicit couch angle list (space-separated)\n"
                      << "  --bixel-width <mm>            Bixel width in mm (default: 7)\n\n"
                      << "Angle semantics:\n"
                      << "  Explicit lists (--gantry-angles + --couch-angles): paired 1:1\n"
                      << "    beam[i] = (gantryAngles[i], couchAngles[i])\n"
                      << "  Start/step/stop with couch-step > 0: Cartesian product (multi-arc)\n"
                      << "    totalBeams = numGantry x numCouch\n"
                      << "    Each couch angle sweeps through ALL gantry angles.\n\n"
                      << "Examples:\n"
                      << "  plan --mode photons --gantry-step 60 --bixel-width 7\n"
                      << "  plan --gantry-start 0 --gantry-stop 360 --gantry-step 90 --couch-start 10\n"
                      << "  plan --gantry-angles 0 90 180 270 --couch-angles 0 10 20 30\n"
                      << "  plan --gantry-step 4 --couch-start -5 --couch-step 5 --couch-stop 10\n"
                      << "    -> 90 gantry x 3 couch = 270 beams\n"
                      << "=============================\n";
        } else if (cmd == "loadPhaseSpace") {
            std::vector<std::string> phspArgs(tokens.begin() + 1, tokens.end());
            loadPhaseSpace(phspArgs, state);
        } else if (cmd == "doseCalc") {
            std::vector<std::string> dcArgs(tokens.begin() + 1, tokens.end());
            doseCalc(dcArgs, state);
        } else if (cmd == "optimize") {
            std::vector<std::string> optArgs(tokens.begin() + 1, tokens.end());
            optimize(optArgs, state);
        } else if (cmd == "leafSeq") {
            std::vector<std::string> lsArgs(tokens.begin() + 1, tokens.end());
            leafSeq(lsArgs, state);
        } else if (cmd == "phsp-info") {
            if (state.phaseSpaceSources.empty()) {
                std::cout << "No phase-space sources loaded.\n";
            } else {
                std::cout << "\n=== Phase-Space Sources (" << state.phaseSpaceSources.size() << " beams) ===\n";
                for (size_t i = 0; i < state.phaseSpaceSources.size(); ++i) {
                    const auto& m = state.phaseSpaceSources[i]->getMetrics();
                    std::cout << "Beam " << i << " (gantry=" << state.phaseSpaceSources[i]->getGantryAngle()
                              << " couch=" << state.phaseSpaceSources[i]->getCouchAngle() << " deg):\n";
                    std::cout << "  Particles: " << m.totalCount
                              << " (photons=" << m.photonCount
                              << " electrons=" << m.electronCount
                              << " positrons=" << m.positronCount << ")\n";
                    std::cout << "  Energy: mean=" << m.meanEnergy << " MeV"
                              << " [" << m.minEnergy << ", " << m.maxEnergy << "]\n";
                    std::cout << "  Angular: sigma_u=" << m.angularSpreadU
                              << " sigma_v=" << m.angularSpreadV << "\n";
                }
                std::cout << "=========================\n";
            }
        } else if (cmd == "info") {
            std::cout << "\n=== Current State ===\n";
            std::cout << "Patient data:   " << (state.patientData ? "loaded" : "not loaded") << "\n";
            std::cout << "Plan:           " << (state.plan ? state.plan->getName() : "not created") << "\n";
            if (state.plan) {
                std::cout << "Machine:        " << state.plan->getMachine().getName()
                          << (state.plan->getMachine().isPhaseSpace() ? " [phase-space]" : " [generic]")
                          << "\n";
            }
            std::cout << "STF:            " << (state.stf ? "generated" : "not generated") << "\n";
            std::cout << "Dij:            " << (state.dij ? "computed (nnz=" + std::to_string(state.dij->getNumNonZeros()) + ")" : "not computed") << "\n";
            std::cout << "Optimization:   " << (!state.optimizedWeights.empty() ? "done" : "not done") << "\n";
            std::cout << "Leaf Sequencing: " << (state.leafSequenceDone() ?
                std::to_string(state.leafSequences.size()) + " beams sequenced" : "not done") << "\n";
            std::cout << "Phase-space:    " << (!state.phaseSpaceSources.empty() ?
                std::to_string(state.phaseSpaceSources.size()) + " beams loaded" : "not loaded") << "\n";
            if (state.stf) {
                state.stf->printSummary();
                std::cout << "\nTo proceed: STF is ready for dose calculation\n";
            } else if (state.plan) {
                std::cout << "\nPlan Details:\n";
                state.plan->printSummary();
                std::cout << "\nTo proceed: use 'loadPhaseSpace' for phase-space machines,\n";
                std::cout << "             or check why STF was not generated\n";
            } else if (state.patientData) {
                std::cout << "\nTo proceed: use 'plan [options]' to generate a treatment plan\n";
                std::cout << "             or 'plan-help' to see available options\n";
            } else {
                std::cout << "\nTo proceed: use 'load <dicom_dir>' to load patient data first\n";
            }
            std::cout << "=====================\n";
        } else {
            std::cerr << "Unknown command: " << cmd << "\n";
        }
    }
    return 0;
}

// ── loadPhaseSpace: builds beam sources from IAEA PSF files ──
int loadPhaseSpace(const std::vector<std::string>& args, WorkflowState& state) {
    if (!state.plan) {
        std::cerr << "Error: No plan generated. Use 'plan [options]' first.\n";
        return 1;
    }

    if (!state.plan->getMachine().isPhaseSpace()) {
        std::cerr << "Error: Current machine is not a phase-space machine.\n";
        std::cerr << "Use 'plan --machine Varian_TrueBeam6MV' to select a PSF machine.\n";
        return 1;
    }

    // Default parameters (collimator shared across beams, couch from plan)
    double collimatorAngle = 0.0;
    int64_t maxParticles = 1000000;
    int64_t vizSample = 100000;

    // Parse options
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--collimator" && i + 1 < args.size()) {
            collimatorAngle = std::stod(args[++i]);
        } else if (args[i] == "--max-particles" && i + 1 < args.size()) {
            maxParticles = std::stoll(args[++i]);
        } else if (args[i] == "--viz-sample" && i + 1 < args.size()) {
            vizSample = std::stoll(args[++i]);
        }
    }

    PhaseSpaceBuildOptions options;
    options.collimatorAngle = collimatorAngle;
    options.maxParticles = maxParticles;
    options.vizSampleSize = vizSample;

    std::cout << "\n=== Loading Phase-Space Beam Sources ===\n";

    try {
        state.phaseSpaceSources = PhaseSpaceBuilder::build(*state.plan, options);

        for (size_t i = 0; i < state.phaseSpaceSources.size(); ++i) {
            const auto& metrics = state.phaseSpaceSources[i]->getMetrics();
            std::cout << "Beam " << i << " (gantry=" << state.phaseSpaceSources[i]->getGantryAngle()
                      << " couch=" << state.phaseSpaceSources[i]->getCouchAngle() << " deg): "
                      << metrics.totalCount << " particles, "
                      << "mean E=" << metrics.meanEnergy << " MeV\n";
        }

        std::cout << "\n=== Loaded " << state.phaseSpaceSources.size() << " beam sources ===\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}

// ── leafSeq command: leaf sequencing from optimized weights ──
int leafSeq(const std::vector<std::string>& args, WorkflowState& state) {
    if (state.optimizedWeights.empty()) {
        std::cerr << "Error: No optimized weights. Use 'optimize' first.\n";
        return 1;
    }

    int numLevels = 15;
    double minSegmentMU = 0.0;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--num-levels" && i + 1 < args.size()) {
            numLevels = std::stoi(args[++i]);
        } else if (args[i] == "--min-segment-mu" && i + 1 < args.size()) {
            minSegmentMU = std::stod(args[++i]);
        }
    }

    std::cout << "\n=== Leaf Sequencing ===\n";
    std::cout << "Levels: " << numLevels << " | Min segment MU: " << minSegmentMU << "\n";

    LeafSequencerOptions opts;
    opts.numLevels = numLevels;
    opts.minSegmentMU = minSegmentMU;

    // Use leaf position resolution from machine geometry if available
    if (state.plan) {
        double res = state.plan->getMachine().getGeometry().leafPositionResolution;
        if (res > 0.0) opts.leafPositionResolution = res;
    }

    auto progressCb = [](int current, int total) {
        std::cout << "\r  Beam " << current << "/" << total << std::flush;
    };

    try {
        auto result = LeafSequencingPipeline::run(
            state.optimizedWeights, *state.stf, *state.dij,
            *state.plan, *state.patientData, *state.doseGrid,
            opts, progressCb);

        state.leafSequences = std::move(result.beamSequences);
        state.deliverableWeights = std::move(result.deliverableWeights);
        state.deliverableStats = std::move(result.deliverableStats);

        std::cout << "\n\nPer-beam summary:\n";
        std::cout << std::setw(6) << "Beam" << std::setw(10) << "Segments"
                  << std::setw(12) << "MU" << std::setw(12) << "Fidelity" << "\n";
        for (size_t i = 0; i < state.leafSequences.size(); ++i) {
            const auto& seq = state.leafSequences[i];
            std::cout << std::setw(6) << i
                      << std::setw(10) << seq.segments.size()
                      << std::setw(12) << std::fixed << std::setprecision(1) << seq.totalMU
                      << std::setw(12) << std::setprecision(4) << seq.fluenceFidelity << "\n";
        }
        std::cout << "\nTotal: " << result.totalSegments << " segments, "
                  << std::fixed << std::setprecision(1) << result.totalMU << " MU, "
                  << "mean fidelity " << std::setprecision(4) << result.meanFidelity << "\n";

        std::cout << "\n--- Deliverable Dose Statistics ---\n";
        PlanAnalysis::print(state.deliverableStats);

    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }

    std::cout << "=== Done ===\n";
    return 0;
}

// ── doseCalc command: compute Dij with optional caching ──
int doseCalc(const std::vector<std::string>& args, WorkflowState& state) {
    if (!state.stf || state.stf->isEmpty()) {
        std::cerr << "Error: No STF generated. Use 'plan' first.\n";
        return 1;
    }

    double doseResolution = 2.5;
    bool useCache = true;
    double absoluteThreshold = 1e-6;
    double relativeThreshold = 0.01;
    int numThreads = 0;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--dose-resolution" && i + 1 < args.size()) {
            doseResolution = std::stod(args[++i]);
        } else if (args[i] == "--no-cache") {
            useCache = false;
        } else if (args[i] == "--abs-threshold" && i + 1 < args.size()) {
            absoluteThreshold = std::stod(args[++i]);
        } else if (args[i] == "--rel-threshold" && i + 1 < args.size()) {
            relativeThreshold = std::stod(args[++i]);
        } else if (args[i] == "--threads" && i + 1 < args.size()) {
            numThreads = std::stoi(args[++i]);
        }
    }

    std::cout << "\n=== Dose Calculation ===\n";
    std::cout << "Resolution: " << doseResolution << " mm | Cache: " << (useCache ? "on" : "off") << "\n";

    DoseCalcPipelineOptions opts;
    opts.resolution = {doseResolution, doseResolution, doseResolution};
    opts.useCache = useCache;
    opts.absoluteThreshold = absoluteThreshold;
    opts.relativeThreshold = relativeThreshold;
    opts.numThreads = numThreads;

    auto progressCb = [](int current, int total) {
        std::cout << "\r  Beam " << current << "/" << total << std::flush;
    };

    try {
        auto result = DoseCalculationPipeline::run(
            *state.plan, *state.stf, *state.patientData, opts, progressCb);
        state.dij = result.dij;
        state.doseGrid = result.doseGrid;
        state.computeGrid = result.doseGrid;

        std::cout << "\nDij: " << state.dij->getNumVoxels() << " x " << state.dij->getNumBixels()
                  << " (nnz: " << state.dij->getNumNonZeros() << ")"
                  << (result.cacheHit ? " [cached]" : "") << "\n";
    } catch (const std::exception& e) {
        std::cerr << "\nError: " << e.what() << "\n";
        return 1;
    }

    std::cout << "=== Done ===\n";
    return 0;
}

// ── optimize command: run L-BFGS-B fluence optimization ──
int optimize(const std::vector<std::string>& args, WorkflowState& state) {
    if (!state.dij || state.dij->getNumNonZeros() == 0) {
        std::cerr << "Error: No Dij computed. Use 'doseCalc' first.\n";
        return 1;
    }

    int maxIter = 400;
    double tolerance = 1e-5;
    double targetDose = 66.0;
    double oarMaxDose = 30.0;
    double smooth = 0.0;
    double l2Reg = 0.0;
    double l1Reg = 0.0;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--max-iter" && i + 1 < args.size()) {
            maxIter = std::stoi(args[++i]);
        } else if (args[i] == "--tolerance" && i + 1 < args.size()) {
            tolerance = std::stod(args[++i]);
        } else if (args[i] == "--target-dose" && i + 1 < args.size()) {
            targetDose = std::stod(args[++i]);
        } else if (args[i] == "--oar-max-dose" && i + 1 < args.size()) {
            oarMaxDose = std::stod(args[++i]);
        } else if (args[i] == "--smooth" && i + 1 < args.size()) {
            smooth = std::stod(args[++i]);
        } else if (args[i] == "--l2-reg" && i + 1 < args.size()) {
            l2Reg = std::stod(args[++i]);
        } else if (args[i] == "--l1-reg" && i + 1 < args.size()) {
            l1Reg = std::stod(args[++i]);
        }
    }

    std::cout << "\n=== Optimization ===\n";
    std::cout << "Max iter: " << maxIter << " | Tolerance: " << tolerance
              << " | Target dose: " << targetDose << " Gy\n";
    if (smooth > 0.0) std::cout << "Spatial smoothing: " << smooth << "\n";
    if (l2Reg > 0.0) std::cout << "L2 regularization: " << l2Reg << "\n";
    if (l1Reg > 0.0) std::cout << "L1 regularization (MU penalty): " << l1Reg << "\n";

    OptimizationConfig config;
    config.maxIterations = maxIter;
    config.tolerance = tolerance;
    config.targetDose = targetDose;
    config.spatialSmoothingWeight = smooth;
    config.l2RegWeight = l2Reg;
    config.l1RegWeight = l1Reg;

    // Store prescribed dose in plan (single source of truth)
    if (state.plan) {
        state.plan->setPrescribedDose(targetDose);
    }

    auto protocol = ObjectiveProtocol::lungIMRT(targetDose);

    try {
        auto result = OptimizationPipeline::run(
            *state.dij, config, protocol, *state.patientData, *state.doseGrid,
            state.stf ? state.stf.get() : nullptr);

        state.optimizedWeights = std::move(result.weights);
        state.doseResult = result.doseResult;

        std::cout << "\nOptimization " << (result.converged ? "converged" : "reached max iterations")
                  << " in " << result.iterations << " iterations\n"
                  << "Final objective: " << result.finalObjective << "\n";

        PlanAnalysis::print(result.stats);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    std::cout << "=== Done ===\n";
    return 0;
}

int main(int argc, char* argv[]) {
    Logger::init();

    WorkflowState state;

    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    std::string command = argv[1];

    if (command == "help" || command == "-h" || command == "--help") {
        printUsage(argv[0]);
        return 0;
    }

    if (command == "load") {
        if (argc < 3) {
            std::cerr << "Error: Missing DICOM directory path\n";
            return 1;
        }
        return loadDicom(argv[2], state);
    }


    if (command == "plan") {
        std::cerr << "Error: Use 'interactive' mode to load data then create a plan.\n";
        std::cerr << "  " << argv[0] << " interactive\n";
        return 1;
    }

    if (command == "generateStf") {
        std::cerr << "Error: 'generateStf' has been merged into 'plan'. STF is now generated automatically.\n";
        return 1;
    }

    if (command == "interactive") {
        return runInteractive(state);
    }

    if (command == "optimize") {
        std::vector<std::string> optArgs;
        for (int i = 2; i < argc; ++i) optArgs.emplace_back(argv[i]);
        return optimize(optArgs, state);
    }

    if (command == "leafSeq") {
        std::vector<std::string> lsArgs;
        for (int i = 2; i < argc; ++i) lsArgs.emplace_back(argv[i]);
        return leafSeq(lsArgs, state);
    }

    if (command == "doseCalc") {
        std::vector<std::string> dcArgs;
        for (int i = 2; i < argc; ++i) dcArgs.emplace_back(argv[i]);
        return doseCalc(dcArgs, state);
    }

    std::cerr << "Unknown command: " << command << "\n";
    printUsage(argv[0]);
    return 1;
}
