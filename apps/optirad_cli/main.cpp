#include "io/DicomImporter.hpp"
#include "core/Patient.hpp"
#include "core/PatientData.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Structure.hpp"
#include "geometry/Volume.hpp"
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
    
    // Full import - loads CT volume, structures with contours, patient info
    auto patientData = importer.importAll(path);
    
    if (!patientData) {
        std::cerr << "Error: Failed to import patient data\n";
        return 1;
    }
    
    // Display patient info
    if (auto* patient = patientData->getPatient()) {
        std::cout << "\nPatient Information:\n";
        std::cout << "  Name: " << patient->getName() << "\n";
        std::cout << "  ID:   " << patient->getID() << "\n";
    }
    
    // Display CT volume info
    if (auto* ct = patientData->getCTVolume()) {
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
        
        // Find HU range
        int16_t minHU = 32767, maxHU = -32768;
        for (size_t i = 0; i < ct->size(); ++i) {
            minHU = std::min(minHU, ct->data()[i]);
            maxHU = std::max(maxHU, ct->data()[i]);
        }
        std::cout << "  HU Range:   [" << minHU << ", " << maxHU << "]\n";
    }
    
    // Display structures
    if (auto* structures = patientData->getStructureSet()) {
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
