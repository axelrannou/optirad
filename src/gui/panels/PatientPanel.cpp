#include "PatientPanel.hpp"

namespace optirad {

void PatientPanel::update() {
    // Update logic
}

void PatientPanel::render() {
    // ImGui rendering
    // if (ImGui::Begin("Patient")) {
    //     if (m_patient) {
    //         ImGui::Text("Name: %s", m_patient->getName().c_str());
    //         ImGui::Text("ID: %s", m_patient->getID().c_str());
    //     }
    //     if (m_structureSet) {
    //         for (size_t i = 0; i < m_structureSet->getCount(); ++i) {
    //             const auto* s = m_structureSet->getStructure(i);
    //             if (s) ImGui::Text("  - %s", s->getName().c_str());
    //         }
    //     }
    // }
    // ImGui::End();
}

} // namespace optirad
