#include "Theme.hpp"

namespace optirad {

// --- Dark theme palette -------------------------------------------------
// Base black + marine (#274678 = RGB 39,70,120) accent
static constexpr ImVec4 kBlack        = {0.00f, 0.00f, 0.00f, 1.00f};
static constexpr ImVec4 kBlackBg      = {0.06f, 0.06f, 0.06f, 1.00f}; // softer black for panels
static constexpr ImVec4 kMarine       = {0.153f, 0.275f, 0.471f, 1.00f}; // #274678
static constexpr ImVec4 kMarineHov    = {0.200f, 0.345f, 0.565f, 1.00f}; // brighter hover
static constexpr ImVec4 kMarineAct    = {0.245f, 0.420f, 0.655f, 1.00f}; // pressed/active
static constexpr ImVec4 kMarineDark   = {0.095f, 0.175f, 0.305f, 1.00f}; // dark variant (docked tab bg)
static constexpr ImVec4 kWhite        = {1.00f, 1.00f, 1.00f, 1.00f};
static constexpr ImVec4 kWhiteDim     = {0.80f, 0.80f, 0.80f, 1.00f};
static constexpr ImVec4 kMarineBorder = {0.18f, 0.32f, 0.53f, 1.00f};

// --- Light theme palette -------------------------------------------------
static constexpr ImVec4 kLightBg      = {1.00f, 1.00f, 1.00f, 1.00f};
static constexpr ImVec4 kLightPane    = {0.96f, 0.96f, 0.96f, 1.00f}; // panel/child bg
static constexpr ImVec4 kLightCtrl    = {0.92f, 0.92f, 0.92f, 1.00f}; // buttons, frames, tabs
static constexpr ImVec4 kLightHov     = {0.85f, 0.85f, 0.85f, 1.00f};
static constexpr ImVec4 kLightAct     = {0.78f, 0.78f, 0.78f, 1.00f};
static constexpr ImVec4 kLightBorder  = {0.75f, 0.75f, 0.75f, 1.00f};
static constexpr ImVec4 kBlackText    = {0.00f, 0.00f, 0.00f, 1.00f};
static constexpr ImVec4 kDarkGrayText = {0.20f, 0.20f, 0.20f, 1.00f};

// --- Specific theme colors -----------------------------------------------

static ThemeColors g_colors;

const ThemeColors& getThemeColors() {
    return g_colors;
}

void applyTheme(AppTheme theme) {
    ImGuiStyle& s = ImGui::GetStyle();
    ImVec4* c = s.Colors;

    s.WindowRounding    = 4.0f;
    s.ChildRounding     = 4.0f;
    s.FrameRounding     = 3.0f;
    s.PopupRounding     = 3.0f;
    s.ScrollbarRounding = 3.0f;
    s.GrabRounding      = 3.0f;
    s.TabRounding       = 3.0f;
    s.FramePadding      = {6.0f, 4.0f};
    s.ItemSpacing       = {8.0f, 5.0f};
    s.ScrollbarSize     = 14.0f;

    if (theme == AppTheme::Dark) {
        // Text
        c[ImGuiCol_Text]                  = kWhite;
        c[ImGuiCol_TextDisabled]          = kWhiteDim;
        c[ImGuiCol_InputTextCursor] = kWhite;

        // Text status
        g_colors.passText = ImVec4(0.4f, 1.0f, 0.4f, 1.0f);
        g_colors.failText = ImVec4(1.0f, 0.2f, 0.2f, 1.0f);
        g_colors.warningText = ImVec4(1.0f, 1.0f, 0.4f, 1.0f);
        g_colors.progressText = ImVec4(1.0f, 1.0f, 0.0f, 1.0f);

        // Generic, phase space specific colours
        g_colors.phaseSpace = ImVec4(0.5f, 0.8f, 1.0f, 1.0f);
        g_colors.generic    = ImVec4(0.8f, 0.8f, 0.5f, 1.0f);

        // DVH plot colors
        g_colors.dvhBackground = ImVec4(0.118f, 0.118f, 0.118f, 1.0f);
        g_colors.dvhGrid       = ImVec4(0.235f, 0.235f, 0.235f, 1.0f);
        g_colors.dvhLabel      = ImVec4(0.784f, 0.784f, 0.784f, 1.0f);
        g_colors.dvhBorder     = ImVec4(0.588f, 0.588f, 0.588f, 1.0f);

        // Windows / containers
        c[ImGuiCol_WindowBg]              = kBlack;
        c[ImGuiCol_ChildBg]               = kBlackBg;
        c[ImGuiCol_PopupBg]               = {0.08f, 0.08f, 0.08f, 0.98f};

        // Borders
        c[ImGuiCol_Border]                = kMarineBorder;
        c[ImGuiCol_BorderShadow]          = kBlack;

        // Frames (InputText, Combo, etc.)
        c[ImGuiCol_FrameBg]               = kMarineDark;
        c[ImGuiCol_FrameBgHovered]        = kMarineHov;
        c[ImGuiCol_FrameBgActive]         = kMarineAct;

        // Title bar
        c[ImGuiCol_TitleBg]               = kMarine;
        c[ImGuiCol_TitleBgActive]         = kMarineAct;
        c[ImGuiCol_TitleBgCollapsed]      = kMarineDark;

        // Menu bar
        c[ImGuiCol_MenuBarBg]             = kMarine;

        // Scrollbar
        c[ImGuiCol_ScrollbarBg]           = kBlackBg;
        c[ImGuiCol_ScrollbarGrab]         = kMarine;
        c[ImGuiCol_ScrollbarGrabHovered]  = kMarineHov;
        c[ImGuiCol_ScrollbarGrabActive]   = kMarineAct;

        // Check mark
        c[ImGuiCol_CheckMark]             = kWhite;

        // Slider / Grab
        c[ImGuiCol_SliderGrab]            = {0.350f, 0.495f, 0.715f, 1.00f};
        c[ImGuiCol_SliderGrabActive]      = {0.450f, 0.595f, 0.815f, 1.00f};

        // Buttons
        c[ImGuiCol_Button]                = kMarine;
        c[ImGuiCol_ButtonHovered]         = kMarineHov;
        c[ImGuiCol_ButtonActive]          = kMarineAct;

        // Headers (selectable, tree node, etc.)
        c[ImGuiCol_Header]                = kMarine;
        c[ImGuiCol_HeaderHovered]         = kMarineHov;
        c[ImGuiCol_HeaderActive]          = kMarineAct;

        // Separator
        c[ImGuiCol_Separator]             = kMarineBorder;
        c[ImGuiCol_SeparatorHovered]      = kMarineHov;
        c[ImGuiCol_SeparatorActive]       = kMarineAct;

        // Resize grip
        c[ImGuiCol_ResizeGrip]            = kMarine;
        c[ImGuiCol_ResizeGripHovered]     = kMarineHov;
        c[ImGuiCol_ResizeGripActive]      = kMarineAct;

        // Tabs
        c[ImGuiCol_Tab]                   = kMarineDark;
        c[ImGuiCol_TabHovered]            = kMarineHov;
        c[ImGuiCol_TabActive]             = kMarine;
        c[ImGuiCol_TabUnfocused]          = kMarineDark;
        c[ImGuiCol_TabUnfocusedActive]    = kMarineDark;

        // Docking
        c[ImGuiCol_DockingPreview]        = kMarineHov;
        c[ImGuiCol_DockingEmptyBg]        = kBlack;

        // Plot
        c[ImGuiCol_PlotLines]             = kWhiteDim;
        c[ImGuiCol_PlotLinesHovered]      = kWhite;
        c[ImGuiCol_PlotHistogram]         = kMarineHov;
        c[ImGuiCol_PlotHistogramHovered]  = kMarineHov;

        // Table
        c[ImGuiCol_TableHeaderBg]         = kMarine;
        c[ImGuiCol_TableBorderStrong]     = kMarineBorder;
        c[ImGuiCol_TableBorderLight]      = kMarineDark;
        c[ImGuiCol_TableRowBg]            = {0.00f, 0.00f, 0.00f, 0.00f};
        c[ImGuiCol_TableRowBgAlt]         = {1.00f, 1.00f, 1.00f, 0.05f};

        // Text selected / drag-drop
        c[ImGuiCol_TextSelectedBg]        = kMarineHov;
        c[ImGuiCol_DragDropTarget]        = kMarineHov;

        // Modal / nav
        c[ImGuiCol_ModalWindowDimBg]      = {0.00f, 0.00f, 0.00f, 0.60f};
        c[ImGuiCol_NavHighlight]          = kMarineHov;
        c[ImGuiCol_NavWindowingHighlight] = kWhiteDim;
        c[ImGuiCol_NavWindowingDimBg]     = {0.00f, 0.00f, 0.00f, 0.40f};

    } else {
        // Text
        c[ImGuiCol_Text]                  = kBlackText;
        c[ImGuiCol_TextDisabled]          = kDarkGrayText;
        c[ImGuiCol_InputTextCursor] = kBlackText;

        // Text status same as dark but more saturated
        g_colors.passText = ImVec4(0.0f, 0.7f, 0.0f, 1.0f);
        g_colors.failText = ImVec4(1.0f, 0.0f, 0.0f, 1.0f);
        g_colors.warningText = ImVec4(0.8f, 0.8f, 0.0f, 1.0f);
        // brighter yellow but more saturated than dark theme for better visibility on light bg
        g_colors.progressText = ImVec4(1.0f, 0.8f, 0.0f, 1.0f);

        // Generic, phase space specific colours
        g_colors.phaseSpace = ImVec4(0.1f, 0.4f, 0.8f, 1.0f);
        g_colors.generic    = ImVec4(0.5f, 0.5f, 0.1f, 1.0f);

        // DVH plot colors
        g_colors.dvhBackground = ImVec4(1.00f, 1.00f, 1.00f, 1.0f);
        g_colors.dvhGrid       = ImVec4(0.85f, 0.85f, 0.85f, 1.0f);
        g_colors.dvhLabel      = ImVec4(0.20f, 0.20f, 0.20f, 1.0f);
        g_colors.dvhBorder     = ImVec4(0.50f, 0.50f, 0.50f, 1.0f);

        // Windows / containers
        c[ImGuiCol_WindowBg]              = kLightBg;
        c[ImGuiCol_ChildBg]               = kLightPane;
        c[ImGuiCol_PopupBg]               = kLightBg;

        // Borders
        c[ImGuiCol_Border]                = kLightBorder;
        c[ImGuiCol_BorderShadow]          = {0.00f, 0.00f, 0.00f, 0.10f};

        // Frames
        c[ImGuiCol_FrameBg]               = kLightCtrl;
        c[ImGuiCol_FrameBgHovered]        = kLightHov;
        c[ImGuiCol_FrameBgActive]         = kLightAct;

        // Title bar
        c[ImGuiCol_TitleBg]               = kLightPane;
        c[ImGuiCol_TitleBgActive]         = kLightCtrl;
        c[ImGuiCol_TitleBgCollapsed]      = kLightBg;

        // Menu bar
        c[ImGuiCol_MenuBarBg]             = kLightCtrl;

        // Scrollbar
        c[ImGuiCol_ScrollbarBg]           = kLightPane;
        c[ImGuiCol_ScrollbarGrab]         = kLightAct;
        c[ImGuiCol_ScrollbarGrabHovered]  = kLightBorder;
        c[ImGuiCol_ScrollbarGrabActive]   = {0.65f, 0.65f, 0.65f, 1.00f};

        // Check mark
        c[ImGuiCol_CheckMark]             = kBlackText;

        // Slider / Grab
        c[ImGuiCol_SliderGrab]            = {0.65f, 0.65f, 0.65f, 1.00f};
        c[ImGuiCol_SliderGrabActive]      = {0.55f, 0.55f, 0.55f, 1.00f};

        // Buttons
        c[ImGuiCol_Button]                = kLightCtrl;
        c[ImGuiCol_ButtonHovered]         = kLightHov;
        c[ImGuiCol_ButtonActive]          = kLightAct;

        // Headers
        c[ImGuiCol_Header]                = kLightCtrl;
        c[ImGuiCol_HeaderHovered]         = kLightHov;
        c[ImGuiCol_HeaderActive]          = kLightAct;

        // Separator
        c[ImGuiCol_Separator]             = kLightBorder;
        c[ImGuiCol_SeparatorHovered]      = kLightAct;
        c[ImGuiCol_SeparatorActive]       = {0.65f, 0.65f, 0.65f, 1.00f};

        // Resize grip
        c[ImGuiCol_ResizeGrip]            = kLightCtrl;
        c[ImGuiCol_ResizeGripHovered]     = kLightHov;
        c[ImGuiCol_ResizeGripActive]      = kLightAct;

        // Tabs
        c[ImGuiCol_Tab]                   = kLightCtrl;
        c[ImGuiCol_TabHovered]            = kLightHov;
        c[ImGuiCol_TabActive]             = kLightBg;
        c[ImGuiCol_TabUnfocused]          = kLightCtrl;
        c[ImGuiCol_TabUnfocusedActive]    = kLightPane;

        // Docking
        c[ImGuiCol_DockingPreview]        = {0.40f, 0.40f, 0.40f, 0.50f};
        c[ImGuiCol_DockingEmptyBg]        = kLightBg;

        // Plot
        c[ImGuiCol_PlotLines]             = kDarkGrayText;
        c[ImGuiCol_PlotLinesHovered]      = kBlackText;
        c[ImGuiCol_PlotHistogram]         = kLightAct;
        c[ImGuiCol_PlotHistogramHovered]  = kLightBorder;

        // Table
        c[ImGuiCol_TableHeaderBg]         = kLightCtrl;
        c[ImGuiCol_TableBorderStrong]     = kLightBorder;
        c[ImGuiCol_TableBorderLight]      = kLightCtrl;
        c[ImGuiCol_TableRowBg]            = {0.00f, 0.00f, 0.00f, 0.00f};
        c[ImGuiCol_TableRowBgAlt]         = {0.00f, 0.00f, 0.00f, 0.04f};

        // Text selected / drag-drop
        c[ImGuiCol_TextSelectedBg]        = {0.65f, 0.65f, 0.65f, 0.40f};
        c[ImGuiCol_DragDropTarget]        = kLightAct;

        // Modal / nav
        c[ImGuiCol_ModalWindowDimBg]      = {0.00f, 0.00f, 0.00f, 0.25f};
        c[ImGuiCol_NavHighlight]          = kLightHov;
        c[ImGuiCol_NavWindowingHighlight] = kLightAct;
        c[ImGuiCol_NavWindowingDimBg]     = {0.00f, 0.00f, 0.00f, 0.20f};
    }
}

} // namespace optirad
