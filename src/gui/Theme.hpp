#pragma once

#include <imgui.h>

namespace optirad {

enum class AppTheme { Dark, Light };

struct ThemeColors {
    ImVec4 phaseSpace;
    ImVec4 generic;
    ImVec4 passText;
    ImVec4 failText;
    ImVec4 warningText;
    ImVec4 progressText;
    
    // DVH plot colors
    ImVec4 dvhBackground;
    ImVec4 dvhGrid;
    ImVec4 dvhLabel;
    ImVec4 dvhBorder;
};

const ThemeColors& getThemeColors();

/// Apply the chosen color theme globally to ImGui::GetStyle().
void applyTheme(AppTheme theme);

} // namespace optirad
