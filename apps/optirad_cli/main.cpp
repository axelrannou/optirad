#include "io/DicomImporter.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "core/Plan.hpp"
#include "core/Machine.hpp"
#include "steering/StfProperties.hpp"
#include "steering/IStfGenerator.hpp"
#include "steering/PhotonIMRTStfGenerator.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <string>
#include <memory>
#include <sstream>
#include <vector>
#include <cstdlib>

using namespace optirad;

// Shared state across commands
struct AppState {
    std::shared_ptr<PatientData> patientData;
    std::shared_ptr<Plan> plan;
};

// Forward declarations
int generateStf(const std::vector<std::string>& args, AppState& state);
std::unique_ptr<optirad::IStfGenerator> selectStfGenerator(const std::string& mode, double gantryStart, double gantryStep, double gantryStop, double bixelWidth, const std::array<double, 3>& iso);

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " <command> [options]\n\n"
              << "Commands:\n"
              << "  load <dicom_dir>                Load and inspect DICOM directory\n"
              << "  plan [options]                   Generate a treatment plan (requires DICOM data first)\n"
              << "  generateStf                      Generate STF properties (requires plan first)\n"
              << "  interactive                      Enter interactive mode\n"
              << "  help                             Show this help message\n\n"
              << "Plan options:\n"
              << "  --mode <photons|protons>         Radiation mode (default: photons)\n"
              << "  --machine <name>                 Machine name (default: Generic)\n"
              << "  --fractions <n>                  Number of fractions (default: 30)\n"
              << "  --gantry-start <deg>             Gantry start angle (default: 0)\n"
              << "  --gantry-step <deg>              Gantry angle step (default: 4)\n"
              << "  --gantry-stop <deg>              Gantry stop angle exclusive (default: 360)\n"
              << "  --bixel-width <mm>               Bixel width (default: 7)\n";
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
        std::cout << "\nGeometric Information:\n";
        std::cout << "  Patient Position: " << grid.getPatientPosition() << "\n";
        std::cout << "  Slice Thickness:  " << grid.getSliceThickness() << " mm\n";
        std::cout << "  Image Orientation (row): ["
                  << orientation[0] << ", " << orientation[1] << ", " << orientation[2] << "]\n";
        std::cout << "  Image Orientation (col): ["
                  << orientation[3] << ", " << orientation[4] << ", " << orientation[5] << "]\n";

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
    double bixelWidth = 7.0;

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

    // Create machine
    if (machineName == "Generic") {
        plan->setMachine(Machine::createGenericPhoton());
    } else {
        std::cerr << "Warning: Unknown machine '" << machineName << "', using Generic.\n";
        plan->setMachine(Machine::createGenericPhoton());
    }

    // Configure STF properties
    StfProperties stf;
    stf.setGantryAngles(gantryStart, gantryStep, gantryStop);
    stf.bixelWidth = bixelWidth;

    // Compute isocenter from Target structures or fallback to CT center
    auto iso = plan->computeIsoCenter();
    stf.setUniformIsoCenter(iso);

    plan->setStfProperties(stf);

    // Store in state
    state.plan = plan;

    // Print summary
    plan->printSummary();

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
                      << "  plan [options]                Create treatment plan (requires DICOM data)\n"
                      << "  plan-help                     Show plan options and examples\n"
                      << "  generateStf                   Generate STF from plan (requires plan first)\n"
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
                      << "  --bixel-width <mm>            Bixel width in mm (default: 7)\n\n"
                      << "Examples:\n"
                      << "  plan --mode photons --gantry-step 60 --bixel-width 7\n"
                      << "  plan --fractions 25 --machine Generic\n"
                      << "  plan --gantry-start 0 --gantry-stop 360 --gantry-step 30\n"
                      << "=============================\n";
        } else if (cmd == "generateStf") {
            std::vector<std::string> stfArgs(tokens.begin() + 1, tokens.end());
            generateStf(stfArgs, state);
        } else if (cmd == "info") {
            std::cout << "\n=== Current State ===\n";
            std::cout << "Patient data: " << (state.patientData ? "loaded" : "not loaded") << "\n";
            std::cout << "Plan:         " << (state.plan ? state.plan->getName() : "not created") << "\n";
            if (state.plan) {
                std::cout << "\nPlan Details:\n";
                state.plan->printSummary();
                std::cout << "\nTo proceed: use 'generateStf' to generate STF from this plan\n";
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

    // Use plan's STF properties to generate STF
    const auto& stfProps = state.plan->getStfProperties();
    std::string radiationMode = state.plan->getRadiationMode();
    
    // Get isocenter from plan
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    if (!stfProps.isoCenters.empty()) {
        iso = stfProps.isoCenters[0];
    }

    // Infer parameters from current STF properties
    double gantryStart = !stfProps.gantryAngles.empty() ? stfProps.gantryAngles[0] : 0.0;
    double gantryStop = !stfProps.gantryAngles.empty() ? stfProps.gantryAngles.back() + 1.0 : 360.0;
    double gantryStep = stfProps.gantryAngles.size() > 1 ? stfProps.gantryAngles[1] - stfProps.gantryAngles[0] : 60.0;
    double bixelWidth = stfProps.bixelWidth;

    auto generator = selectStfGenerator(radiationMode, gantryStart, gantryStep, gantryStop, bixelWidth, iso);
    auto stf = generator->generate();

    // Print summary
    std::cout << "\n=== Generated STF Properties (from plan) ===\n";
    std::cout << "  Radiation Mode: " << radiationMode << "\n";
    std::cout << "  Gantry angles: ";
    for (auto a : stf->gantryAngles) std::cout << a << " ";
    std::cout << "\n  Couch angles: ";
    for (auto a : stf->couchAngles) std::cout << a << " ";
    std::cout << "\n  Bixel width: " << stf->bixelWidth << " mm\n";
    std::cout << "  Number of beams: " << stf->numOfBeams << "\n";
    if (!stf->isoCenters.empty()) {
        std::cout << "  IsoCenters: [";
        for (const auto& iso : stf->isoCenters) {
            std::cout << "(" << iso[0] << "," << iso[1] << "," << iso[2] << ") ";
        }
        std::cout << "]\n";
    }
    std::cout << "=========================================\n";
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
        std::vector<std::string> stfArgs;
        for (int i = 2; i < argc; ++i) stfArgs.emplace_back(argv[i]);
        return generateStf(stfArgs, state);
    }

    if (command == "interactive") {
        return runInteractive(state);
    }

    if (command == "optimize") {
        std::cout << "Optimization not yet implemented\n";
        return 0;
    }

    std::cerr << "Unknown command: " << command << "\n";
    printUsage(argv[0]);
    return 1;
}
