#pragma once

#include "Patient.hpp"
#include "ImportedPlan.hpp"
#include "geometry/Volume.hpp"
#include "geometry/StructureSet.hpp"
#include "geometry/Grid.hpp"
#include <array>
#include <memory>
#include <optional>
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
    
    // Plan imported from a DICOM RTPLAN (nullptr when the patient has none)
    void setImportedPlan(std::shared_ptr<const ImportedPlan> plan) { m_importedPlan = std::move(plan); }
    const ImportedPlan* getImportedPlan() const { return m_importedPlan.get(); }
    bool hasImportedPlan() const { return m_importedPlan != nullptr && !m_importedPlan->beams.empty(); }

    // Voxel spacing (x, y, z mm) of the imported RTDOSE grid, used as the default dose-grid resolution
    void setImportedDoseGridResolution(const std::array<double, 3>& res) { m_importedDoseResolution = res; }
    const std::optional<std::array<double, 3>>& getImportedDoseGridResolution() const { return m_importedDoseResolution; }

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
    std::shared_ptr<const ImportedPlan> m_importedPlan;
    std::optional<std::array<double, 3>> m_importedDoseResolution;
    Grid m_emptyGrid;
};

} // namespace optirad
