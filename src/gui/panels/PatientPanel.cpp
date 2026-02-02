#include "PatientPanel.hpp"

namespace optirad {

void PatientPanel::update() {
    // Update logic
}

void PatientPanel::render() {
    if (!m_visible) return;
    
    // TODO: ImGui::Begin("Patient");
    // ImGui::Text("Patient ID: %s", m_patient.getId().c_str());
    // ImGui::Text("Patient Name: %s", m_patient.getName().c_str());
    // ImGui::Separator();
    // ImGui::Text("Structures: %zu", m_structureSet.getNumStructures());
    // for (size_t i = 0; i < m_structureSet.getNumStructures(); ++i) {
    //     ImGui::Text("  - %s", m_structureSet.getStructure(i).getName().c_str());
    // }
    // ImGui::End();
}

void PatientPanel::setPatient(const Patient& patient) { m_patient = patient; }
void PatientPanel::setStructureSet(const StructureSet& structureSet) { m_structureSet = structureSet; }

} // namespace optirad
