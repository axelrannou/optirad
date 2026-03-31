#include "Application.hpp"
#include "Theme.hpp"
#include "../../external/stb_image.h"
#include "panels/PatientPanel.hpp"
#include "panels/PlanningPanel.hpp"
#include "panels/StfPanel.hpp"
#include "panels/PhaseSpacePanel.hpp"
#include "panels/OptimizationPanel.hpp"
#include "panels/DoseStatsPanel.hpp"
#include "panels/DVHPanel.hpp"
#include "views/SliceView.hpp"
#include "views/View3D.hpp"
#include "views/renderers/PhaseSpaceRenderer.hpp"
#include "utils/Logger.hpp"

#include <imgui.h>
#include <imgui_internal.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>

namespace optirad {

Application::Application() = default;
Application::~Application() = default;

bool Application::loadTexture(const std::string& path, GLuint* outTexture) {
    int w, h, channels;
    stbi_set_flip_vertically_on_load(false);
    unsigned char* pixels = stbi_load(path.c_str(), &w, &h, &channels, 4);
    if (!pixels) {
        Logger::error("loadTexture: failed to load '" + path + "': " +
                      std::string(stbi_failure_reason()));
        return false;
    }
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    glBindTexture(GL_TEXTURE_2D, 0);
    stbi_image_free(pixels);
    *outTexture = tex;
    return true;
}

bool Application::init() {
    // Initialize GLFW
    if (!glfwInit()) {
        Logger::error("Failed to initialize GLFW");
        return false;
    }
    
    // Create window
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA

    // Detect primary monitor resolution (supports 4K etc.)
    const GLFWvidmode* mode = glfwGetVideoMode(glfwGetPrimaryMonitor());
    int initW = mode ? mode->width  : 1920;
    int initH = mode ? mode->height : 1080;
    m_window = glfwCreateWindow(initW, initH, "OptiRad TPS", nullptr, nullptr);
    if (!m_window) {
        Logger::error("Failed to create GLFW window");
        glfwTerminate();
        return false;
    }
    
    glfwMaximizeWindow(m_window);
    glfwMakeContextCurrent(m_window);
    glfwSwapInterval(1); // VSync

    // Initialize GLEW
    glewExperimental = GL_TRUE;
    if (glewInit() != GLEW_OK) {
        Logger::error("Failed to initialize GLEW");
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }

    glEnable(GL_MULTISAMPLE); // Enable MSAA
    
    // Initialize ImGui
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    
    ImGuiIO& io = ImGui::GetIO();
    io.ConfigFlags |= ImGuiConfigFlags_DockingEnable;
    
    // Apply initial dark theme
    applyTheme(m_theme);

    // DPI scaling
    // On Linux/X11, glfwGetWindowContentScale often returns the desktop scale
    // factor (e.g. 2.0 on 4K) but the framebuffer is already at native resolution,
    // so we query framebuffer vs window size to detect actual content scaling.
    float xscale = 1.0f, yscale = 1.0f;
    glfwGetWindowContentScale(m_window, &xscale, &yscale);
    int fbW, fbH, winW, winH;
    glfwGetFramebufferSize(m_window, &fbW, &fbH);
    glfwGetWindowSize(m_window, &winW, &winH);
    // If framebuffer == window size, the OS is NOT doing pixel-doubling,
    // so content scale > 1 likely just means a high-DPI panel at native res.
    if (winW > 0 && fbW > 0) {
        float actualScale = static_cast<float>(fbW) / static_cast<float>(winW);
        if (actualScale < 1.5f) {
            // No pixel doubling — use 1.0 so we don't double-scale
            xscale = 1.0f;
        }
    }
    m_dpiScale = xscale;
    
    // Use a clean, readable font size
    float fontSize = 20.0f * m_dpiScale;
    ImFontConfig fontConfig;
    fontConfig.SizePixels = fontSize;
    io.Fonts->AddFontDefault(&fontConfig);
    
    if (!ImGui_ImplGlfw_InitForOpenGL(m_window, true)) {
        Logger::error("Failed to initialize ImGui GLFW backend");
        ImGui::DestroyContext();
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }
    
    if (!ImGui_ImplOpenGL3_Init("#version 330")) {
        Logger::error("Failed to initialize ImGui OpenGL3 backend");
        ImGui_ImplGlfw_Shutdown();
        ImGui::DestroyContext();
        glfwDestroyWindow(m_window);
        glfwTerminate();
        return false;
    }
    
    // Load theme toggle icon
    loadTexture(OPTIRAD_GUI_IMG_DIR "/day-and-night.png", &m_themeIconTexture);

    // Create 3D view
    m_view3D = std::make_unique<View3D>();
    m_view3D->init();
    m_view3D->setDarkMode(m_theme == AppTheme::Dark);
    
    // Create panels
    m_patientPanel = std::make_unique<PatientPanel>(m_appState);
    m_planningPanel = std::make_unique<PlanningPanel>(m_appState);
    m_stfPanel = std::make_unique<StfPanel>(m_appState);
    m_phaseSpacePanel = std::make_unique<PhaseSpacePanel>(m_appState);
    m_optimizationPanel = std::make_unique<OptimizationPanel>(m_appState);
    m_doseStatsPanel = std::make_unique<DoseStatsPanel>(m_appState);
    m_dvhPanel = std::make_unique<DVHPanel>(m_appState);
    
    // Phase Space panel hidden by default
    m_phaseSpacePanel->setVisible(false);
    
    // Create slice views
    m_axialView = std::make_unique<SliceView>(SliceOrientation::Axial);
    m_sagittalView = std::make_unique<SliceView>(SliceOrientation::Sagittal);
    m_coronalView = std::make_unique<SliceView>(SliceOrientation::Coronal);
    
    return true;
}

void Application::setupDockLayout() {
    ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
    
    // Always clear and rebuild the layout when this is called
    ImGui::DockBuilderRemoveNode(dockspace_id);
    ImGui::DockBuilderAddNode(dockspace_id, ImGuiDockNodeFlags_DockSpace);
    
    ImGuiViewport* viewport = ImGui::GetMainViewport();
    ImGui::DockBuilderSetNodeSize(dockspace_id, viewport->WorkSize);
    ImGui::DockBuilderSetNodePos(dockspace_id, viewport->WorkPos);
    
    // Split: left 35% for sidebars | right 65% for views + bottom
    ImGuiID leftNode, rightNode;
    ImGui::DockBuilderSplitNode(dockspace_id, ImGuiDir_Left, 0.35f, &leftNode, &rightNode);
    
    // Split left into two sidebars: sidebar1 (50% of left = 15%) | sidebar2 (50% of left = 15%)
    ImGuiID sidebar1Node, sidebar2Node;
    ImGui::DockBuilderSplitNode(leftNode, ImGuiDir_Left, 0.50f, &sidebar1Node, &sidebar2Node);
    
    // Split right into views (top 75%) and bottom (25% for dose stats/DVH)
    ImGuiID viewsNode, bottomNode;
    ImGui::DockBuilderSplitNode(rightNode, ImGuiDir_Down, 0.25f, &bottomNode, &viewsNode);
    
    // Split views area into 2x2 grid
    // First split horizontal: top row | bottom row
    ImGuiID topRow, bottomRow;
    ImGui::DockBuilderSplitNode(viewsNode, ImGuiDir_Down, 0.50f, &bottomRow, &topRow);
    
    // Split top row: Axial (left) | Sagittal (right)
    ImGuiID axialNode, sagittalNode;
    ImGui::DockBuilderSplitNode(topRow, ImGuiDir_Left, 0.50f, &axialNode, &sagittalNode);
    
    // Split bottom row: 3D View (left, under Axial) | Coronal (right)
    ImGuiID view3dNode, coronalNode;
    ImGui::DockBuilderSplitNode(bottomRow, ImGuiDir_Left, 0.50f, &view3dNode, &coronalNode);
    
    // Split sidebar1 vertically: Patient (top 35%) | Planning (bottom 65%)
    ImGuiID sidebar1Top, sidebar1Bottom;
    ImGui::DockBuilderSplitNode(sidebar1Node, ImGuiDir_Down, 0.65f, &sidebar1Bottom, &sidebar1Top);
    
    // Split sidebar2 vertically: STF (top 35%) | Optimization (bottom 65%)
    ImGuiID sidebar2Top, sidebar2Bottom;
    ImGui::DockBuilderSplitNode(sidebar2Node, ImGuiDir_Down, 0.65f, &sidebar2Bottom, &sidebar2Top);
    
    // Dock windows into nodes
    // Sidebar 1: Patient (top) + Planning (bottom)
    ImGui::DockBuilderDockWindow("Patient Data", sidebar1Top);
    ImGui::DockBuilderDockWindow("Planning", sidebar1Bottom);
    
    // Sidebar 2: STF (top) + Optimization (bottom)
    ImGui::DockBuilderDockWindow("STF Generation", sidebar2Top);
    ImGui::DockBuilderDockWindow("Optimization", sidebar2Bottom);
    
    // Phase Space (hidden but docked in sidebar2 STF node as a tab)
    ImGui::DockBuilderDockWindow("Phase Space", sidebar2Top);
    
    // Views
    ImGui::DockBuilderDockWindow("Axial", axialNode);
    ImGui::DockBuilderDockWindow("Sagittal", sagittalNode);
    ImGui::DockBuilderDockWindow("Coronal", coronalNode);
    ImGui::DockBuilderDockWindow("3D View", view3dNode);
    ImGui::DockBuilderDockWindow("DVH", view3dNode);  // tabbed with 3D View
    
    // Bottom: Dose Statistics
    ImGui::DockBuilderDockWindow("Dose Statistics", bottomNode);
    
    ImGui::DockBuilderFinish(dockspace_id);
}

void Application::run() {
    while (!glfwWindowShouldClose(m_window)) {
        glfwPollEvents();
        
        // Clear default framebuffer
        int display_w, display_h;
        glfwGetFramebufferSize(m_window, &display_w, &display_h);
        glViewport(0, 0, display_w, display_h);
        glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Update data bridges
        if (m_patientPanel->hasData()) {
            auto* data = m_patientPanel->getPatientData();

            if (!m_appState.dicomLoaded()) {
                m_appState.patientData = std::shared_ptr<PatientData>(
                    data, [](PatientData*) {});

                // Add imported RT Dose to DoseManager (once)
                if (!m_importedDoseAdded && data->hasImportedDose()) {
                    m_appState.doseManager.addDose(
                        "Imported (DICOM)",
                        data->getImportedDose(),
                        data->getImportedDoseGrid());
                    m_appState.syncSelectedDose();
                    m_importedDoseAdded = true;
                }
            }

            m_view3D->setPatientData(data);
            m_axialView->setPatientData(data);
            m_sagittalView->setPatientData(data);
            m_coronalView->setPatientData(data);
        }

        if (m_appState.stfGenerated()) {
            m_view3D->setStf(m_appState.stf.get());
            m_stfPanel->setBeamRenderer(m_view3D->getBeamRenderer());
        }

        if (m_appState.phaseSpaceLoaded()) {
            std::vector<const PhaseSpaceBeamSource*> sources;
            sources.reserve(m_appState.phaseSpaceSources.size());
            for (const auto& src : m_appState.phaseSpaceSources) {
                sources.push_back(src.get());
            }
            m_view3D->setPhaseSpaceSources(sources);
            m_phaseSpacePanel->setPhaseSpaceRenderer(m_view3D->getPhaseSpaceRenderer());
        }
        
        // Start ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // Re-apply theme every frame (ImGui immediate-mode best practice)
        applyTheme(m_theme);
        
        // Create fullscreen dockspace host window
        ImGuiViewport* viewport = ImGui::GetMainViewport();
        ImGui::SetNextWindowPos(viewport->WorkPos);
        ImGui::SetNextWindowSize(viewport->WorkSize);
        ImGui::SetNextWindowViewport(viewport->ID);
        
        ImGuiWindowFlags hostFlags = 
            ImGuiWindowFlags_NoDocking | ImGuiWindowFlags_NoTitleBar |
            ImGuiWindowFlags_NoCollapse | ImGuiWindowFlags_NoResize |
            ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoBringToFrontOnFocus |
            ImGuiWindowFlags_NoNavFocus | ImGuiWindowFlags_NoBackground |
            ImGuiWindowFlags_MenuBar;
        
        ImGui::PushStyleVar(ImGuiStyleVar_WindowRounding, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0.0f);
        ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0.0f, 0.0f));
        
        ImGui::Begin("DockSpaceHost", nullptr, hostFlags);
        ImGui::PopStyleVar(3);
        
        // Render menu bar inside host window
        renderMenuBar();
        
        // Setup initial layout BEFORE the first DockSpace call
        ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
        if (!m_layoutInitialized) {
            setupDockLayout();
            m_layoutInitialized = true;
        }
        
        // Create the dockspace
        ImGui::DockSpace(dockspace_id, ImVec2(0.0f, 0.0f), ImGuiDockNodeFlags_None);
        
        ImGui::End(); // DockSpaceHost
        
        // Render panels
        m_patientPanel->render();
        m_planningPanel->render();
        m_stfPanel->render();
        if (m_phaseSpacePanel->isVisible()) {
            m_phaseSpacePanel->render();
        }
        m_optimizationPanel->render();
        m_doseStatsPanel->render();
        
        // Pass dose data to slice views AFTER panels render
        if (m_appState.doseAvailable() && m_appState.doseResult && m_appState.doseGrid) {
            m_axialView->setDoseData(m_appState.doseResult.get(), m_appState.doseGrid.get());
            m_sagittalView->setDoseData(m_appState.doseResult.get(), m_appState.doseGrid.get());
            m_coronalView->setDoseData(m_appState.doseResult.get(), m_appState.doseGrid.get());
        } else {
            m_axialView->setDoseData(nullptr, nullptr);
            m_sagittalView->setDoseData(nullptr, nullptr);
            m_coronalView->setDoseData(nullptr, nullptr);
        }
        
        // Render slice views (each creates its own ImGui window)
        m_axialView->render();
        m_sagittalView->render();
        m_coronalView->render();
        
        // Render 3D View inside an ImGui window
        render3DViewWindow();
        m_dvhPanel->render();

        //Auto-switch to 3D View for first render
        if (!focus3DView) {
            ImGui::SetWindowFocus("3D View");
            focus3DView = true;
        }

        // Auto-switch to DVH tab when optimization finishes
        if (m_appState.optimizationJustFinished.exchange(false)) {
            ImGui::SetWindowFocus("DVH");
        }

        // Finalize ImGui frame
        ImGui::Render();
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        
        glfwSwapBuffers(m_window);
    }
}

void Application::render3DViewWindow() {
    ImGui::PushStyleVar(ImGuiStyleVar_WindowPadding, ImVec2(0, 0));
    if (!m_view3DVisible) {
        ImGui::PopStyleVar();
        return;
    }
    ImGui::Begin("3D View", &m_view3DVisible);
    ImGui::PopStyleVar();
    
    ImVec2 avail = ImGui::GetContentRegionAvail();
    int w = static_cast<int>(avail.x);
    int h = static_cast<int>(avail.y);
    
    if (w > 0 && h > 0) {
        // Render the 3D scene to FBO
        m_view3D->render(w, h);
        
        // Display the FBO texture (flipped vertically since OpenGL coordinates are bottom-up)
        ImGui::Image(
            reinterpret_cast<void*>(static_cast<intptr_t>(m_view3D->getTextureID())),
            avail,
            ImVec2(0, 1), ImVec2(1, 0) // UV flipped for OpenGL
        );
        
        // Handle mouse input when 3D view is hovered
        if (ImGui::IsItemHovered()) {
            m_view3D->handleImGuiInput();
        }
    }
    
    ImGui::End();
}

void Application::shutdown() {
    if (m_view3D) {
        m_view3D->cleanup();
    }
    
    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    
    if (m_themeIconTexture) {
        glDeleteTextures(1, &m_themeIconTexture);
        m_themeIconTexture = 0;
    }
    if (m_window) {
        glfwDestroyWindow(m_window);
    }
    glfwTerminate();
}

void Application::renderMenuBar() {
    if (ImGui::BeginMenuBar()) {
        if (ImGui::BeginMenu("File")) {
            if (ImGui::MenuItem("Import DICOM...")) {
                // Trigger import dialog in patient panel
            }
            ImGui::Separator();
            if (ImGui::MenuItem("Exit")) {
                glfwSetWindowShouldClose(m_window, true);
            }
            ImGui::EndMenu();
        }
        
        if (ImGui::BeginMenu("View")) {
            // Toggle panel visibility
            bool patientVis = m_patientPanel->isVisible();
            if (ImGui::MenuItem("Patient Panel", nullptr, &patientVis)) {
                m_patientPanel->setVisible(patientVis);
            }
            
            bool planningVis = m_planningPanel->isVisible();
            if (ImGui::MenuItem("Planning Panel", nullptr, &planningVis)) {
                m_planningPanel->setVisible(planningVis);
            }
            
            bool stfVis = m_stfPanel->isVisible();
            if (ImGui::MenuItem("STF Panel", nullptr, &stfVis)) {
                m_stfPanel->setVisible(stfVis);
            }
            
            bool phaseVis = m_phaseSpacePanel->isVisible();
            if (ImGui::MenuItem("Phase Space Panel", nullptr, &phaseVis)) {
                m_phaseSpacePanel->setVisible(phaseVis);
            }
            
            bool optVis = m_optimizationPanel->isVisible();
            if (ImGui::MenuItem("Optimization Panel", nullptr, &optVis)) {
                m_optimizationPanel->setVisible(optVis);
            }
            
            bool doseVis = m_doseStatsPanel->isVisible();
            if (ImGui::MenuItem("Dose Statistics", nullptr, &doseVis)) {
                m_doseStatsPanel->setVisible(doseVis);
            }

            bool dvhVis = m_dvhPanel->isVisible();
            if (ImGui::MenuItem("DVH", nullptr, &dvhVis)) {
                m_dvhPanel->setVisible(dvhVis);
            }
            
            ImGui::Separator();
            
            bool axialVis = m_axialView->isVisible();
            if (ImGui::MenuItem("Axial View", nullptr, &axialVis)) {
                m_axialView->setVisible(axialVis);
            }
            
            bool sagVis = m_sagittalView->isVisible();
            if (ImGui::MenuItem("Sagittal View", nullptr, &sagVis)) {
                m_sagittalView->setVisible(sagVis);
            }
            
            bool corVis = m_coronalView->isVisible();
            if (ImGui::MenuItem("Coronal View", nullptr, &corVis)) {
                m_coronalView->setVisible(corVis);
            }
            
            bool view3DVis = m_view3DVisible;
            if (ImGui::MenuItem("3D View", nullptr, &view3DVis)) {
                m_view3DVisible = view3DVis;
            }
            
            ImGui::Separator();
            
            if (ImGui::MenuItem("Reset Layout")) {
                m_layoutInitialized = false;
                //make all panels visible again
                m_patientPanel->setVisible(true);
                m_planningPanel->setVisible(true);
                m_stfPanel->setVisible(true);
                m_phaseSpacePanel->setVisible(false);
                m_optimizationPanel->setVisible(true);
                m_doseStatsPanel->setVisible(true);
                m_dvhPanel->setVisible(true);
                m_axialView->setVisible(true);
                m_sagittalView->setVisible(true);
                m_coronalView->setVisible(true);
                m_view3DVisible = true;
                // Force rebuild on next frame by removing the existing node
                ImGuiID dockspace_id = ImGui::GetID("MainDockSpace");
                ImGui::DockBuilderRemoveNode(dockspace_id);
            }
            
            ImGui::EndMenu();
        }

        // --- Theme toggle (after separator, right side of menu bar) ---
        ImGui::Separator();

        // Draw the icon button
        ImVec2 iconSize(20.0f * m_dpiScale, 20.0f * m_dpiScale);
        ImGui::PushStyleColor(ImGuiCol_Button,        ImVec4(0, 0, 0, 0));
        ImGui::PushStyleColor(ImGuiCol_ButtonHovered, ImVec4(0.3f, 0.3f, 0.3f, 0.4f));
        ImGui::PushStyleColor(ImGuiCol_ButtonActive,  ImVec4(0.5f, 0.5f, 0.5f, 0.4f));
        if (m_themeIconTexture) {
            if (ImGui::ImageButton("##themeToggle",
                    reinterpret_cast<void*>(static_cast<intptr_t>(m_themeIconTexture)),
                    iconSize)) {
                m_theme = (m_theme == AppTheme::Dark) ? AppTheme::Light : AppTheme::Dark;
                applyTheme(m_theme);
                m_view3D->setDarkMode(m_theme == AppTheme::Dark);
            }
        } else {
            // Fallback text button if icon didn't load
            if (ImGui::Button(m_theme == AppTheme::Dark ? "Light" : "Dark")) {
                m_theme = (m_theme == AppTheme::Dark) ? AppTheme::Light : AppTheme::Dark;
                applyTheme(m_theme);
                m_view3D->setDarkMode(m_theme == AppTheme::Dark);
            }
        }
        ImGui::PopStyleColor(3);
        if (ImGui::IsItemHovered()) {
            ImGui::SetTooltip("%s mode", m_theme == AppTheme::Dark ? "Switch to Light" : "Switch to Dark");
        }

        ImGui::EndMenuBar();
    }
}

} // namespace optirad
