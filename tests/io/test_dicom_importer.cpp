#include "io/DicomImporter.hpp"
#include "core/Patient.hpp"
#include "geometry/StructureSet.hpp"
#include "utils/Logger.hpp"

#include <iostream>
#include <filesystem>

using namespace optirad;

int main(int argc, char* argv[]) {
    std::cout << "=== OptiRad DICOM Importer Test ===\n\n";
    
    Logger::init();
    
    std::string testPath = ".";
    if (argc > 1) {
        testPath = argv[1];
    } else {
        std::cout << "Usage: " << argv[0] << " <dicom_directory>\n";
        std::cout << "No path provided, using current directory.\n\n";
    }
    
    DicomImporter importer;
    
    // Test canImport
    std::cout << "Testing canImport(\"" << testPath << "\")...\n";
    bool canImport = importer.canImport(testPath);
    std::cout << "  Result: " << (canImport ? "YES" : "NO") << "\n\n";
    
    if (!canImport) {
        std::cout << "Cannot import from this path.\n";
        std::cout << "Make sure the directory contains .dcm files.\n";
        return 1;
    }
    
    // Test loadDirectory
    std::cout << "Loading DICOM directory...\n";
    bool loaded = importer.loadDirectory(testPath);
    std::cout << "  Result: " << (loaded ? "SUCCESS" : "FAILED") << "\n\n";
    
    if (!loaded) {
        std::cout << "Failed to load DICOM directory.\n";
        return 1;
    }
    
    // Test importPatient
    std::cout << "Importing patient information...\n";
    auto patient = importer.importPatient(testPath);
    if (patient) {
        std::cout << "  Patient Name: " << patient->getName() << "\n";
        std::cout << "  Patient ID:   " << patient->getID() << "\n\n";
    } else {
        std::cout << "  No patient data found.\n\n";
    }
    
    // Test importStructures
    std::cout << "Importing structures...\n";
    auto structures = importer.importStructures(testPath);
    if (structures) {
        std::cout << "  Found " << structures->getCount() << " structures.\n";
    } else {
        std::cout << "  No structures found.\n";
    }
    
    std::cout << "\n=== Test Complete ===\n";
    return 0;
}
