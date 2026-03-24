#pragma once

#include "Patient.hpp"
#include "geometry/Volume.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Grid.hpp"
#include "dose/DoseMatrix.hpp"
#include <memory>
#include <vector>

namespace optirad {

/**
 * Central container for all patient-related data
 * Analogous to matRad's ct + cst structures
 */
class PatientData {
public:
    PatientData() = default;
    
    // Patient info
    void setPatient(std::unique_ptr<Patient> patient) { m_patient = std::move(patient); }
    Patient* getPatient() { return m_patient.get(); }
    const Patient* getPatient() const { return m_patient.get(); }
    
    // CT volume (HU values stored as int16_t)
    void setCTVolume(std::unique_ptr<Volume<int16_t>> ct) { m_ctVolume = std::move(ct); }
    Volume<int16_t>* getCTVolume() { return m_ctVolume.get(); }
    const Volume<int16_t>* getCTVolume() const { return m_ctVolume.get(); }
    
    // CT in electron density (for dose calculation)
    void setEDVolume(std::unique_ptr<Volume<double>> ed) { m_edVolume = std::move(ed); }
    Volume<double>* getEDVolume() { return m_edVolume.get(); }
    const Volume<double>* getEDVolume() const { return m_edVolume.get(); }
    
    // Structures
    void setStructureSet(std::unique_ptr<StructureSet> structures) { m_structures = std::move(structures); }
    StructureSet* getStructureSet() { return m_structures.get(); }
    const StructureSet* getStructureSet() const { return m_structures.get(); }
    
    // Imported RT Dose (from DICOM)
    void setImportedDose(std::shared_ptr<DoseMatrix> dose, std::shared_ptr<Grid> grid) {
        m_importedDose = std::move(dose);
        m_importedDoseGrid = std::move(grid);
    }
    std::shared_ptr<DoseMatrix> getImportedDose() const { return m_importedDose; }
    std::shared_ptr<Grid> getImportedDoseGrid() const { return m_importedDoseGrid; }
    bool hasImportedDose() const { return m_importedDose != nullptr; }

    // Grid information
    const Grid& getGrid() const { return m_ctVolume ? m_ctVolume->getGrid() : m_emptyGrid; }
    
    // Validation
    bool hasValidCT() const { return m_ctVolume != nullptr && m_ctVolume->size() > 0; }
    bool hasStructures() const { return m_structures != nullptr && m_structures->getCount() > 0; }
    bool isValid() const { return hasValidCT(); }
    
    // HU to ED conversion (simple linear for now)
    void convertHUtoED();
    
private:
    std::unique_ptr<Patient> m_patient;
    std::unique_ptr<Volume<int16_t>> m_ctVolume;   // Hounsfield Units
    std::unique_ptr<Volume<double>> m_edVolume;    // Electron Density
    std::unique_ptr<StructureSet> m_structures;
    std::shared_ptr<DoseMatrix> m_importedDose;
    std::shared_ptr<Grid> m_importedDoseGrid;
    Grid m_emptyGrid;
};

} // namespace optirad
