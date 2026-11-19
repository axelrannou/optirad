#include "DeliveryPanel.hpp"
#include <imgui.h>

namespace optirad {

void DeliveryPanel::renderControls(DeliveryPlayback& delivery) {
    if (!delivery.available()) {
        ImGui::TextWrapped("%s", delivery.description().c_str());
        return;
    }
    if (delivery.playing) {
        const float size = ImGui::GetFrameHeight();
        if (ImGui::Button("##playArc", ImVec2(size, size))) delivery.togglePlaying();
        const ImVec2 min = ImGui::GetItemRectMin();
        const float barW = size * 0.18f, gap = size * 0.16f, top = size * 0.25f, bottom = size * 0.75f;
        const float left = min.x + (size - 2.0f * barW - gap) * 0.5f;
        auto* draw = ImGui::GetWindowDrawList();
        const ImU32 color = ImGui::GetColorU32(ImGuiCol_Text);
        draw->AddRectFilled({left, min.y + top}, {left + barW, min.y + bottom}, color);
        draw->AddRectFilled({left + barW + gap, min.y + top}, {left + 2.0f * barW + gap, min.y + bottom}, color);
    } else if (ImGui::ArrowButton("##playArc", ImGuiDir_Right)) {
        delivery.togglePlaying();
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip(delivery.playing ? "Pause delivery" : "Play whole plan");
    ImGui::SameLine();
    if (ImGui::ArrowButton("##rewindArc", ImGuiDir_Left)) delivery.rewind();
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Rewind whole plan");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(130.0f);
    ImGui::SliderFloat("Playback", &delivery.playbackSpeed, 0.1f, 10.0f, "%.1fx");
    double time = delivery.time();
    const double min = 0.0, max = delivery.duration();
    if (ImGui::SliderScalar("Time", ImGuiDataType_Double, &time, &min, &max, "%.2f s"))
        delivery.seek(time);
    int group = static_cast<int>(delivery.frame().group);
    if (ImGui::SliderInt(delivery.isVmat() ? "Arc" : "Beam", &group, 0,
                         static_cast<int>(delivery.groupCount()) - 1))
        delivery.selectGroup(static_cast<size_t>(group));
    const auto& frame = delivery.frame();
    ImGui::Text("CP %zu | Gantry %.1f | Couch %.1f | %.2f MU/s | %.1f mm/s",
                frame.controlPoint, frame.beam.getGantryAngle(), frame.beam.getCouchAngle(),
                frame.doseRate, frame.leafSpeed);
    ImGui::Text("Gantry speed: %.2f deg/s (%s)", frame.gantrySpeed,
                delivery.estimated() ? "estimated" : "delivery");
    if (!delivery.isVmat())
        ImGui::Text("%s", frame.beamOn ? "Beam on: shooting" : "Beam off: repositioning");
    else if (!frame.beamOn)
        ImGui::TextDisabled("Beam off: repositioning");
}

void DeliveryPanel::render() {
    if (!m_visible) return;
    ImGui::Begin("Delivery", &m_visible);
    m_state.delivery.synchronize(m_state);
    ImGui::Checkbox("Show delivery overlays", &m_state.delivery.showOverlays);
    ImGui::TextWrapped("Yellow overlays in 3D and slice views. Playback is shared with BEV.");
    ImGui::Separator();
    renderControls(m_state.delivery);
    if (m_state.delivery.available()) {
        ImGui::Separator();
        ImGui::TextWrapped("%s", m_state.delivery.description().c_str());
        ImGui::TextDisabled("%s", m_state.delivery.frame().hasAperture ?
            "MLC aperture at source-collimator distance" : "Field outline only (no MLC results)");
    }
    ImGui::End();
}

} // namespace optirad
