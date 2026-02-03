#include "BeamPanel.hpp"

namespace optirad {

void BeamPanel::update() {}

void BeamPanel::render() {
    if (!m_visible) return;
    
    // TODO: ImGui::Begin("Beams");
    // if (m_plan) {
    //     for (size_t i = 0; i < m_plan->getNumBeams(); ++i) {
    //         const auto& beam = m_plan->getBeams()[i];
    //         ImGui::Text("Beam %zu: Gantry=%.1f°, Couch=%.1f°", 
    //                     i, beam.getGantryAngle(), beam.getCouchAngle());
    //     }
    //     if (ImGui::Button("Add Beam")) { /* ... */ }
    // }
    // ImGui::End();
}

void BeamPanel::setPlan(Plan* plan) { m_plan = plan; }

} // namespace optirad
