#pragma once

#include <imgui.h>
#include <algorithm>
#include <cstddef>
#include <cstdio>

namespace optirad::layout {

// Keep columns content-sized; scroll instead of squeezing or restoring stale widths.
inline constexpr ImGuiTableFlags contentTableFlags =
    ImGuiTableFlags_SizingFixedFit | ImGuiTableFlags_ScrollX |
    ImGuiTableFlags_NoSavedSettings;

inline float framedTextWidth(const char* text) {
    return ImGui::CalcTextSize(text, nullptr, true).x + ImGui::GetStyle().FramePadding.x * 2.0f;
}

inline float comboWidth(const char* text) {
    return framedTextWidth(text) + ImGui::GetFrameHeight();
}

inline bool combo(const char* label, int* currentItem, const char* const items[], int count) {
    float width = 0.0f;
    for (int i = 0; i < count; ++i)
        width = std::max(width, comboWidth(items[i]));
    ImGui::SetNextItemWidth(width);
    return ImGui::Combo(label, currentItem, items, count);
}

inline float floatInputWidth(float value) {
    char text[64];
    std::snprintf(text, sizeof(text), "%.1f", value);
    return framedTextWidth(text);
}

inline ImVec2 tableSize(std::size_t rows, float rowHeight, float headerHeight = 0.0f) {
    const auto& style = ImGui::GetStyle();
    return ImVec2(0, std::max(headerHeight, ImGui::GetTextLineHeight()) + style.CellPadding.y * 2.0f +
                     static_cast<float>(rows) * (rowHeight + style.CellPadding.y * 2.0f) +
                     style.ScrollbarSize + 2.0f);
}

inline void sameLineIfFits(float nextWidth) {
    const float nextX = ImGui::GetItemRectMax().x + ImGui::GetStyle().ItemSpacing.x;
    const float right = ImGui::GetWindowPos().x + ImGui::GetWindowContentRegionMax().x;
    if (nextX + nextWidth <= right)
        ImGui::SameLine();
}

inline ImVec2 fullWidthButtonSize(const char* label, float height = 0.0f) {
    return ImVec2(std::max(ImGui::GetContentRegionAvail().x, framedTextWidth(label)),
                  height > 0.0f ? std::max(height, ImGui::GetFrameHeight()) : 0.0f);
}

} // namespace optirad::layout
