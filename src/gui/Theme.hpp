#pragma once

#include <imgui.h>

namespace optirad {

enum class AppTheme { Dark, Light };

/// Apply the chosen color theme globally to ImGui::GetStyle().
void applyTheme(AppTheme theme);

} // namespace optirad
