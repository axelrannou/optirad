#include "io/DicomImporter.hpp"
#include "io/MachineLoader.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Machine.hpp"
#include "core/Stf.hpp"
#include "steering/StfProperties.hpp"
#include "steering/IStfGenerator.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"
#include "phsp/PhaseSpaceBeamSource.hpp"
#include "phsp/PhaseSpaceData.hpp"
#include "dose/DoseEngineFactory.hpp"
#include "dose/DoseInfluenceMatrix.hpp"
#include "dose/DoseMatrix.hpp"
#include "dose/DijSerializer.hpp"
#include "dose/DoseCalcOptions.hpp"
#include "optimization/OptimizerFactory.hpp"
#include "optimization/optimizers/LBFGSOptimizer.hpp"
#include "optimization/objectives/SquaredDeviation.hpp"
#include "optimization/objectives/SquaredOverdose.hpp"
#include "optimization/objectives/SquaredUnderdose.hpp"
#include "optimization/objectives/DVHObjective.hpp"
#include "optimization/Constraint.hpp"
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

using namespace optirad;

// Shared state across commands
struct AppState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
    std::shared_ptr<StfProperties> stfProps;  // lightweight STF properties
    std::shared_ptr<Stf> stf;                 // full STF with beams and rays
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> phaseSpaceSources;

    // Dose calculation
    std::shared_ptr<DoseInfluenceMatrix> dij;
    std::shared_ptr<Grid> doseGrid;

    // Optimization
    std::vector<double> optimizedWeights;
    std::shared_ptr<DoseMatrix> doseResult;
};

// Forward declarations
int loadPhaseSpace(const std::vector<std::string>& args, AppState& state);
int doseCalc(const std::vector<std::string>& args, AppState& state);
int optimize(const std::vector<std::string>& args, AppState& state);
std::unique_ptr<optirad::IStfGenerator> selectStfGenerator(const std::string& mode, double gantryStart, double gantryStep, double gantryStop, double bixelWidth, const std::array<double, 3>& iso);

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " <command> [options]\n\n"
              << "Commands:\n"
              << "  load <dicom_dir>                Load and inspect DICOM directory\n"
              << "  plan [options]                   Generate a treatment plan and STF (requires DICOM data first)\n"
              << "  doseCalc [options]               Calculate Dij (requires plan)\n"
              << "  optimize [options]               Run optimization (requires Dij)\n"
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
              << "Phase-space options:\n"
              << "  --collimator <deg>               Collimator angle (default: 0)\n"
              << "  --max-particles <n>              Max particles per beam (default: 1000000)\n"
              << "  --viz-sample <n>                 Visualization sample per beam (default: 100000)\n"
              << "  (Couch angles come from plan's --couch-start/step/stop or --couch-angles)\n";
}

int loadDicom(const std::string& path, AppState& state) {
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
    state.plan.reset();
    state.stfProps.reset();
    state.stf.reset();
    state.dij.reset();
    state.doseGrid.reset();
    state.optimizedWeights.clear();
    state.doseResult.reset();
    state.phaseSpaceSources.clear();

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

int createPlan(const std::vector<std::string>& args, AppState& state) {
    if (!state.patientData) {
        std::cerr << "Error: No patient data loaded. Use 'load <dicom_dir>' first.\n";
        return 1;
    }

    // Default plan parameters
    std::string radiationMode = "photons";
    std::string machineName = "Generic";
    int numFractions = 30;
    double gantryStart = 0.0;
    double gantryStep = 4.0;
    double gantryStop = 360.0;
    double couchStart = 0.0;
    double couchStep = 0.0;
    double couchStop = 0.0;
    double bixelWidth = 7.0;
    std::vector<double> gantryAnglesList;
    std::vector<double> couchAnglesList;

    // Parse plan options
    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--mode" && i + 1 < args.size()) {
            radiationMode = args[++i];
        } else if (args[i] == "--machine" && i + 1 < args.size()) {
            machineName = args[++i];
        } else if (args[i] == "--fractions" && i + 1 < args.size()) {
            numFractions = std::stoi(args[++i]);
        } else if (args[i] == "--gantry-start" && i + 1 < args.size()) {
            gantryStart = std::stod(args[++i]);
        } else if (args[i] == "--gantry-step" && i + 1 < args.size()) {
            gantryStep = std::stod(args[++i]);
        } else if (args[i] == "--gantry-stop" && i + 1 < args.size()) {
            gantryStop = std::stod(args[++i]);
        } else if (args[i] == "--couch-start" && i + 1 < args.size()) {
            couchStart = std::stod(args[++i]);
        } else if (args[i] == "--couch-step" && i + 1 < args.size()) {
            couchStep = std::stod(args[++i]);
        } else if (args[i] == "--couch-stop" && i + 1 < args.size()) {
            couchStop = std::stod(args[++i]);
        } else if (args[i] == "--gantry-angles") {
            // Consume subsequent numeric arguments
            while (i + 1 < args.size() && args[i + 1][0] != '-') {
                gantryAnglesList.push_back(std::stod(args[++i]));
            }
        } else if (args[i] == "--couch-angles") {
            // Consume subsequent numeric arguments
            while (i + 1 < args.size() && args[i + 1][0] != '-') {
                couchAnglesList.push_back(std::stod(args[++i]));
            }
        } else if (args[i] == "--bixel-width" && i + 1 < args.size()) {
            bixelWidth = std::stod(args[++i]);
        }
    }

    // Create the plan
    auto plan = std::make_shared<Plan>();
    plan->setName("TreatmentPlan");
    plan->setRadiationMode(radiationMode);
    plan->setNumOfFractions(numFractions);
    plan->setPatientData(state.patientData);

    // Load machine from JSON data file
    try {
        plan->setMachine(MachineLoader::load(radiationMode, machineName));
    } catch (const std::exception& e) {
        std::cerr << "Error loading machine '" << machineName << "': " << e.what() << "\n";
        return 1;
    }

    // Configure STF properties
    StfProperties stf;
    if (!gantryAnglesList.empty()) {
        stf.setGantryAngles(gantryAnglesList);
    } else {
        stf.setGantryAngles(gantryStart, gantryStep, gantryStop);
    }
    if (!couchAnglesList.empty()) {
        stf.setCouchAngles(couchAnglesList);
    } else if (couchStep > 0.0) {
        stf.setCouchAngles(couchStart, couchStep, couchStop);
    } else {
        stf.setUniformCouchAngle(couchStart);
    }
    stf.ensureConsistentAngles();
    stf.bixelWidth = bixelWidth;

    // Compute isocenter from Target structures
    auto iso = plan->computeIsoCenter();
    stf.setUniformIsoCenter(iso);

    plan->setStfProperties(stf);

    // Store in state
    state.plan = plan;

    // New plan invalidates previously generated steering/dose artifacts
    state.stfProps.reset();
    state.stf.reset();
    state.dij.reset();
    state.doseGrid.reset();
    state.optimizedWeights.clear();
    state.doseResult.reset();
    state.phaseSpaceSources.clear();

    // Print summary
    plan->printSummary();

    // Automatically generate STF for generic (non-phase-space) machines
    if (!plan->getMachine().isPhaseSpace()) {
        std::cout << "\nGenerating STF...\n";

        const auto& stfPropsRef = plan->getStfProperties();
        std::array<double, 3> isoRef = {0.0, 0.0, 0.0};
        if (!stfPropsRef.isoCenters.empty()) {
            isoRef = stfPropsRef.isoCenters[0];
        }

        PhotonIMRTStfGenerator generator(0.0, 360.0, 360.0, stfPropsRef.bixelWidth, isoRef);
        generator.setMachine(plan->getMachine());
        generator.setRadiationMode(plan->getRadiationMode());
        generator.setGantryAngles(stfPropsRef.gantryAngles);
        generator.setCouchAngles(stfPropsRef.couchAngles);

        auto patientData = plan->getPatientData();
        if (patientData && patientData->hasValidCT() && patientData->hasStructures()) {
            const auto* ct = patientData->getCTVolume();
            const auto* structureSet = patientData->getStructureSet();
            const auto& grid = ct->getGrid();

            generator.setGrid(grid);
            generator.setStructureSet(*structureSet);
            generator.setCTResolution(grid.getSpacing());
        }

        state.stfProps = generator.generate();
        state.stf = std::make_shared<Stf>(generator.generateStf());

        state.stf->printSummary();
    } else {
        std::cout << "\nPhase-space machine: use 'loadPhaseSpace' to build beam sources.\n";
    }

    return 0;
}

int runInteractive(AppState& state) {
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

// Helper: select STF generator based on mode (future extensibility)
std::unique_ptr<optirad::IStfGenerator> selectStfGenerator(const std::string& mode, double gantryStart, double gantryStep, double gantryStop, double bixelWidth, const std::array<double, 3>& iso) {
    // For now, only Photon IMRT is implemented
    if (mode == "photons" || mode == "photon_imrt") {
        return std::make_unique<optirad::PhotonIMRTStfGenerator>(gantryStart, gantryStep, gantryStop, bixelWidth, iso);
    }
    // Future: add more modalities here
    // Default fallback
    return std::make_unique<optirad::PhotonIMRTStfGenerator>(gantryStart, gantryStep, gantryStop, bixelWidth, iso);
}

// New command: generateStf (uses plan parameters)
int generateStf(const std::vector<std::string>& args, AppState& state) {
    if (!state.plan) {
        std::cerr << "Error: No plan generated. Use 'plan [options]' to create a plan first.\n";
        return 1;
    }

    // Gate: STF generation only for generic machines
    if (state.plan->getMachine().isPhaseSpace()) {
        std::cerr << "Error: STF generation is not applicable for phase-space machines.\n";
        std::cerr << "Use 'loadPhaseSpace' instead to build the beam source from IAEA PSF files.\n";
        return 1;
    }

    // Use plan's STF properties to generate STF
    const auto& stfProps = state.plan->getStfProperties();
    std::string radiationMode = state.plan->getRadiationMode();
    
    // Get isocenter from plan
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    if (!stfProps.isoCenters.empty()) {
        iso = stfProps.isoCenters[0];
    }

    double bixelWidth = stfProps.bixelWidth;

    // Create generator and pass the already-expanded paired gantry/couch lists.
    // The plan's StfProperties already performed Cartesian product expansion,
    // so we pass the full paired lists directly to avoid re-inferring start/step/stop.
    PhotonIMRTStfGenerator generator(0.0, 360.0, 360.0, bixelWidth, iso);
    generator.setMachine(state.plan->getMachine());
    generator.setRadiationMode(radiationMode);
    generator.setGantryAngles(stfProps.gantryAngles);
    generator.setCouchAngles(stfProps.couchAngles);

    // Extract target voxel world coordinates from patient data
    auto patientData = state.plan->getPatientData();
    if (patientData && patientData->hasValidCT() && patientData->hasStructures()) {
        const auto* ct = patientData->getCTVolume();
        const auto* structureSet = patientData->getStructureSet();
        const auto& grid = ct->getGrid();

        // Pass Grid and StructureSet to the generator for proper 3D margin expansion
        generator.setGrid(grid);
        generator.setStructureSet(*structureSet);

        // Pass CT resolution for potential padding
        auto spacing = grid.getSpacing();
        generator.setCTResolution(spacing);

        std::cout << "Using Grid and StructureSet for 3D margin expansion\n";
    } else {
        std::cout << "Warning: No patient data available, using fixed field size\n";
    }

    // Also store lightweight StfProperties for backward compatibility
    state.stfProps = generator.generate();

    // Generate full Stf with beams and rays
    auto stf = std::make_shared<Stf>(generator.generateStf());
    state.stf = stf;

    // Print summary
    state.stf->printSummary();
    return 0;
}

// New command: loadPhaseSpace (builds beam sources from IAEA PSF files for all gantry angles)
int loadPhaseSpace(const std::vector<std::string>& args, AppState& state) {
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

    // Get gantry/couch angles and isocenter from plan
    const auto& stfProps = state.plan->getStfProperties();
    const auto& gantryAngles = stfProps.gantryAngles;
    const auto& couchAngles = stfProps.couchAngles;
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    if (!stfProps.isoCenters.empty()) {
        iso = stfProps.isoCenters[0];
    }

    std::cout << "\n=== Loading Phase-Space Beam Sources ===\n";
    std::cout << "Beams: " << gantryAngles.size() << " gantry/couch angle pairs\n";
    std::cout << "Collimator: " << collimatorAngle << " deg\n";
    std::cout << "Max particles/beam: " << maxParticles << "\n";
    std::cout << "Viz sample/beam: " << vizSample << "\n\n";

    try {
        const int numBeams = static_cast<int>(gantryAngles.size());
        state.phaseSpaceSources.clear();
        state.phaseSpaceSources.resize(numBeams);

        // Build all beams in parallel (OpenMP)
        #pragma omp parallel for schedule(dynamic)
        for (int i = 0; i < numBeams; ++i) {
            auto source = std::make_shared<PhaseSpaceBeamSource>();
            source->configure(state.plan->getMachine(), gantryAngles[i],
                              collimatorAngle,
                              (static_cast<size_t>(i) < couchAngles.size() ? couchAngles[i] : 0.0),
                              iso);
            source->build(maxParticles, vizSample);
            state.phaseSpaceSources[i] = std::move(source);
        }

        // Print metrics sequentially
        for (int i = 0; i < numBeams; ++i) {
            const auto& metrics = state.phaseSpaceSources[i]->getMetrics();
            double beamCouch = (static_cast<size_t>(i) < couchAngles.size() ? couchAngles[i] : 0.0);
            std::cout << "Beam " << i << " (gantry=" << gantryAngles[i]
                      << " couch=" << beamCouch << " deg): "
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

// ── doseCalc command: compute Dij with optional caching ──
int doseCalc(const std::vector<std::string>& args, AppState& state) {
    if (!state.stf || state.stf->isEmpty()) {
        std::cerr << "Error: No STF generated. Use 'generateStf' first.\n";
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
    std::cout << "Dose grid resolution: " << doseResolution << " mm\n";
    std::cout << "Thresholds: absolute=" << absoluteThreshold
              << ", relative=" << (relativeThreshold * 100.0) << "%\n";
    std::cout << "Threads: " << (numThreads > 0 ? std::to_string(numThreads) : "all") << "\n";

    // Set dose resolution on plan
    state.plan->setDoseGridResolution({doseResolution, doseResolution, doseResolution});

    // Create dose grid
    const auto& ctGrid = state.patientData->getGrid();
    auto doseGrid = std::make_shared<Grid>(
        Grid::createDoseGrid(ctGrid, {doseResolution, doseResolution, doseResolution}));
    state.doseGrid = doseGrid;

    auto dims = doseGrid->getDimensions();
    std::cout << "Dose grid: " << dims[0] << "x" << dims[1] << "x" << dims[2]
              << " (" << doseGrid->getNumVoxels() << " voxels)\n";

    // Check cache
    if (useCache) {
        std::string patientName = "unknown";
        if (state.patientData->getPatient()) {
            patientName = state.patientData->getPatient()->getName();
        }
        std::string cacheFile = DijSerializer::getCacheDir() + "/" +
            DijSerializer::buildCacheKey(
                patientName,
                static_cast<int>(state.stf->getCount()),
                state.plan->getStfProperties().bixelWidth,
                doseResolution);

        if (DijSerializer::exists(cacheFile)) {
            std::cout << "Loading Dij from cache: " << cacheFile << "\n";
            auto loaded = DijSerializer::load(cacheFile);
            state.dij = std::make_shared<DoseInfluenceMatrix>(std::move(loaded));
            std::cout << "Dij loaded: " << state.dij->getNumVoxels() << " x "
                      << state.dij->getNumBixels() << " (nnz: " << state.dij->getNumNonZeros() << ")\n";
            std::cout << "=== Done ===\n";
            return 0;
        }
    }

    // Compute Dij
    std::cout << "Computing Dij (" << state.stf->getCount() << " beams, "
              << state.stf->getTotalNumOfBixels() << " bixels)...\n";

    auto engine = DoseEngineFactory::create("PencilBeam");
    engine->setProgressCallback([](int current, int total, const std::string& msg) {
        std::cout << "\r  Beam " << current << "/" << total << " " << msg << std::flush;
    });

    // Apply user options
    DoseCalcOptions opts;
    opts.absoluteThreshold = absoluteThreshold;
    opts.relativeThreshold = relativeThreshold;
    opts.numThreads = numThreads;
    engine->setOptions(opts);

    auto start = std::chrono::steady_clock::now();
    auto dij = engine->calculateDij(*state.plan, *state.stf, *state.patientData, *doseGrid);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    state.dij = std::make_shared<DoseInfluenceMatrix>(std::move(dij));
    std::cout << "\nDij: " << state.dij->getNumVoxels() << " x " << state.dij->getNumBixels()
              << " (nnz: " << state.dij->getNumNonZeros() << ") in " << elapsed << "s\n";

    // Save to cache
    if (useCache) {
        std::string patientName = "unknown";
        if (state.patientData->getPatient()) {
            patientName = state.patientData->getPatient()->getName();
        }
        std::string cacheFile = DijSerializer::getCacheDir() + "/" +
            DijSerializer::buildCacheKey(
                patientName,
                static_cast<int>(state.stf->getCount()),
                state.plan->getStfProperties().bixelWidth,
                doseResolution);
        DijSerializer::save(*state.dij, cacheFile);
        std::cout << "Saved to cache: " << cacheFile << "\n";
    }

    std::cout << "=== Done ===\n";
    return 0;
}

// ── optimize command: run L-BFGS-B fluence optimization ──
int optimize(const std::vector<std::string>& args, AppState& state) {
    if (!state.dij || state.dij->getNumNonZeros() == 0) {
        std::cerr << "Error: No Dij computed. Use 'doseCalc' first.\n";
        return 1;
    }

    int maxIter = 400;
    double tolerance = 1e-5;
    double targetDose = 66.0;
    double oarMaxDose = 30.0;

    for (size_t i = 0; i < args.size(); ++i) {
        if (args[i] == "--max-iter" && i + 1 < args.size()) {
            maxIter = std::stoi(args[++i]);
        } else if (args[i] == "--tolerance" && i + 1 < args.size()) {
            tolerance = std::stod(args[++i]);
        } else if (args[i] == "--target-dose" && i + 1 < args.size()) {
            targetDose = std::stod(args[++i]);
        } else if (args[i] == "--oar-max-dose" && i + 1 < args.size()) {
            oarMaxDose = std::stod(args[++i]);
        }
    }

    // Build objectives from structures
    const auto* ss = state.patientData->getStructureSet();
    if (!ss || ss->getCount() == 0) {
        std::cerr << "Error: No structures available for optimization.\n";
        return 1;
    }

    const auto& ctGrid = state.patientData->getGrid();
    const auto& doseGrid = *state.doseGrid;

    // ── Pre-optimization summary ──
    std::cout << "\n╔══════════════════════════════════════════════════════════════════════════════╗\n";
    std::cout << "║                       PRE-OPTIMIZATION SUMMARY                              ║\n";
    std::cout << "╚══════════════════════════════════════════════════════════════════════════════╝\n\n";
    std::cout << "  Bixels:      " << state.dij->getNumBixels() << "\n";
    std::cout << "  Dose voxels: " << state.dij->getNumVoxels() << "\n";
    std::cout << "  NNZ:         " << state.dij->getNumNonZeros() << "\n";
    std::cout << "  Max iter:    " << maxIter << "\n";
    std::cout << "  Tolerance:   " << tolerance << "\n\n";

    std::vector<std::unique_ptr<ObjectiveFunction>> ownedObjs;
    std::vector<ObjectiveFunction*> objPtrs;

    std::cout << "  Objectives:\n";
    std::cout << "  ──────────────────────────────────────────────────────────────\n";

    // Helper lambda to find a structure by name (exact) or pattern (substring)
    auto findStructures = [&](const std::string& pattern, bool exact) 
        -> std::vector<std::pair<const Structure*, std::vector<size_t>>> {
        std::vector<std::pair<const Structure*, std::vector<size_t>>> results;
        for (size_t si = 0; si < ss->getCount(); ++si) {
            const auto* structure = ss->getStructure(si);
            if (!structure || structure->getVoxelIndices().empty()) continue;
            std::string name = structure->getName();
            bool match = exact ? (name == pattern) : (name.find(pattern) != std::string::npos);
            if (match) {
                auto doseIndices = Grid::mapVoxelIndices(ctGrid, doseGrid, structure->getVoxelIndices());
                results.push_back({structure, std::move(doseIndices)});
            }
        }
        return results;
    };

    // Helper lambda to add an objective
    auto addObjective = [&](std::unique_ptr<ObjectiveFunction> obj, 
                            const Structure* structure,
                            const std::vector<size_t>& doseIndices,
                            const std::string& desc) {
        obj->setStructure(structure);
        obj->setVoxelIndices(doseIndices);
        std::cout << "  " << desc << " (" << doseIndices.size() << " voxels)\n";
        objPtrs.push_back(obj.get());
        ownedObjs.push_back(std::move(obj));
    };

    // ── PTVs: MinDVH D98% at targetDose (only PTV N and PTVT) ──
    // PTV ok and PTV_66dosi are unions of PTV N + PTVT, so we only optimize the individual PTVs.
    for (const std::string& ptvName : {"PTV N", "PTVT"}) {
        for (auto& [structure, doseIndices] : findStructures(ptvName, true)) {
            auto obj = std::make_unique<DVHObjective>();
            obj->setType(DVHObjective::Type::MIN_DVH);
            obj->setDoseThreshold(targetDose);
            obj->setVolumeFraction(0.98);
            obj->setWeight(100.0);
            addObjective(std::move(obj), structure, doseIndices,
                "Target: " + ptvName + " -> MinDVH D98% >= " + std::to_string(targetDose) + " Gy (w=100)");
        }
    }

    // ── Poumon_D: SquaredOverdose at MLD 20 Gy ──
    for (auto& [structure, doseIndices] : findStructures("Poumon_D", true)) {
        auto obj = std::make_unique<SquaredOverdose>();
        obj->setMaxDose(20.0);
        obj->setWeight(10.0);
        addObjective(std::move(obj), structure, doseIndices,
            "OAR:    Poumon_D -> SquaredOverdose @ 20 Gy (w=10)");
    }

    // ── Poumon_G: SquaredOverdose at MLD 20 Gy ──
    for (auto& [structure, doseIndices] : findStructures("Poumon_G", true)) {
        auto obj = std::make_unique<SquaredOverdose>();
        obj->setMaxDose(20.0);
        obj->setWeight(10.0);
        addObjective(std::move(obj), structure, doseIndices,
            "OAR:    Poumon_G -> SquaredOverdose @ 20 Gy (w=10)");
    }

    // ── Coeur: MaxDVH V40 < 30% ──
    for (auto& [structure, doseIndices] : findStructures("Coeur", true)) {
        auto obj = std::make_unique<DVHObjective>();
        obj->setType(DVHObjective::Type::MAX_DVH);
        obj->setDoseThreshold(40.0);
        obj->setVolumeFraction(0.30);
        obj->setWeight(30.0);
        addObjective(std::move(obj), structure, doseIndices,
            "OAR:    Coeur -> MaxDVH V40Gy <= 30% (w=30)");
    }

    // ── Oesophage: MaxDVH V50 < 30% ──
    for (auto& [structure, doseIndices] : findStructures("Oesophage", true)) {
        auto obj = std::make_unique<DVHObjective>();
        obj->setType(DVHObjective::Type::MAX_DVH);
        obj->setDoseThreshold(50.0);
        obj->setVolumeFraction(0.30);
        obj->setWeight(50.0);
        addObjective(std::move(obj), structure, doseIndices,
            "OAR:    Oesophage -> MaxDVH V50Gy <= 30% (w=50)");
    }

    // ── Canal_Medullaire: MaxDVH V10Gy <= 0% ──
    for (auto& [structure, doseIndices] : findStructures("Canal_Medullaire", true)) {
        auto obj = std::make_unique<DVHObjective>();
        obj->setType(DVHObjective::Type::MAX_DVH);
        obj->setDoseThreshold(10.0);
        obj->setVolumeFraction(0.0);
        obj->setWeight(200.0);
        addObjective(std::move(obj), structure, doseIndices,
            "OAR:    Canal_Medullaire -> MaxDVH V10Gy <= 0% (w=200)");
    }

    if (objPtrs.empty()) {
        std::cerr << "Error: No valid objectives generated from structures.\n";
        return 1;
    }

    std::cout << "\n  Total objectives: " << objPtrs.size() << "\n";
    std::cout << "\nRunning L-BFGS-B optimizer...\n";

    auto optimizer = OptimizerFactory::create("LBFGS");
    optimizer->setMaxIterations(maxIter);
    optimizer->setTolerance(tolerance);

    // Configure LBFGS-specific parameters
    if (auto* lbfgs = dynamic_cast<LBFGSOptimizer*>(optimizer.get())) {
        // maxFluence must be large enough for the optimizer to reach the prescribed dose.
        // Estimate: targetDose / maxBixelDose. Use a generous margin.
        double maxBixelDose = state.dij->getMaxValue();
        double estimatedMaxFluence = (maxBixelDose > 1e-10)
            ? targetDose / maxBixelDose * 2.0  // 2x headroom
            : 1000.0;
        lbfgs->setMaxFluence(estimatedMaxFluence);
        lbfgs->setPrescriptionDose(targetDose);
        lbfgs->setHotspotThreshold(1.04);   // penalise > 104% of Rx
        lbfgs->setHotspotPenalty(2000.0);

        std::cout << "\n  LBFGS parameters:\n";
        std::cout << "    Max fluence:      " << estimatedMaxFluence
                  << " (estimated from maxBixelDose=" << maxBixelDose << ")\n";
        std::cout << "    Prescription:     " << targetDose << " Gy\n";
        std::cout << "    Hotspot limit:    " << targetDose * 1.04 << " Gy (104%)\n";
        std::cout << "    Hotspot penalty:  2000\n";
    } else {
        std::cerr << "Warning: dynamic_cast to LBFGSOptimizer failed!\n";
    }

    auto start = std::chrono::steady_clock::now();
    std::vector<Constraint> constraints;
    auto result = optimizer->optimize(*state.dij, objPtrs, constraints);
    auto end = std::chrono::steady_clock::now();
    double elapsed = std::chrono::duration<double>(end - start).count();

    state.optimizedWeights = std::move(result.weights);

    std::cout << "\nOptimization " << (result.converged ? "converged" : "reached max iterations")
              << " in " << result.iterations << " iterations (" << elapsed << "s)\n"
              << "Final objective: " << result.finalObjective << "\n";

    // Forward dose computation
    std::cout << "\nComputing forward dose...\n";
    auto engine = DoseEngineFactory::create("PencilBeam");
    auto dose = engine->calculateDose(*state.dij, state.optimizedWeights, *state.doseGrid);
    state.doseResult = std::make_shared<DoseMatrix>(std::move(dose));

    // Comprehensive Plan Analysis
    auto stats = PlanAnalysis::computeStats(
        *state.doseResult, *state.patientData, *state.doseGrid, targetDose);
    PlanAnalysis::print(stats);

    std::cout << "=== Done ===\n";
    return 0;
}

int main(int argc, char* argv[]) {
    Logger::init();

    AppState state;

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

    if (command == "doseCalc") {
        std::vector<std::string> dcArgs;
        for (int i = 2; i < argc; ++i) dcArgs.emplace_back(argv[i]);
        return doseCalc(dcArgs, state);
    }

    std::cerr << "Unknown command: " << command << "\n";
    printUsage(argv[0]);
    return 1;
}
