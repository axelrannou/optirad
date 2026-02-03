#include "io/DicomImporter.hpp"
#include "core/Patient.hpp"
#include "geometry/StructureSet.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <string>

using namespace optirad;

void printUsage(const char* progName) {
    std::cout << "Usage: " << progName << " <command> [options]\n\n"
              << "Commands:\n"
              << "  load <dicom_dir>    Load and inspect DICOM directory\n"
              << "  optimize <config>   Run optimization from config file\n"
              << "  help                Show this help message\n";
}

int loadDicom(const std::string& path) {
    std::cout << "=== OptiRad DICOM Loader ===\n\n";
    std::cout << "Loading DICOM from: " << path << "\n\n";
    
    DicomImporter importer;
    
    if (!importer.canImport(path)) {
        std::cerr << "Error: Cannot import from path: " << path << "\n";
        return 1;
    }
    
    // Load directory to scan for files
    if (!importer.loadDirectory(path)) {
        std::cerr << "Error: Failed to load DICOM directory\n";
        return 1;
    }
    
    // Import patient info
    auto patient = importer.importPatient(path);
    if (patient) {
        std::cout << "Patient Information:\n";
        std::cout << "  Name: " << patient->getName() << "\n";
        std::cout << "  ID:   " << patient->getID() << "\n\n";
    }
    
    // Import structures
    auto structures = importer.importStructures(path);
    if (structures) {
        std::cout << "Structures found: " << structures->getCount() << "\n";
        for (size_t i = 0; i < structures->getCount(); ++i) {
            const auto* s = structures->getStructure(i);
            if (s) {
                std::cout << "  - " << s->getName() << " (" << s->getType() << ")\n";
            }
        }
    }
    
    std::cout << "\nDICOM load complete.\n";
    return 0;
}

int main(int argc, char* argv[]) {
    Logger::init();
    
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
        return loadDicom(argv[2]);
    }
    
    if (command == "optimize") {
        std::cout << "Optimization not yet implemented\n";
        return 0;
    }
    
    std::cerr << "Unknown command: " << command << "\n";
    printUsage(argv[0]);
    return 1;
}
