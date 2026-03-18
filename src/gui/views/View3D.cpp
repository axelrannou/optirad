#include "View3D.hpp"
#include "renderers/CubeRenderer.hpp"
#include "renderers/AxisLabels.hpp"
#include "renderers/VolumeRenderer.hpp"
#include "renderers/StructureRenderer.hpp"
#include "renderers/BeamRenderer.hpp"
#include "renderers/PhaseSpaceRenderer.hpp"
#include "core/PatientData.hpp"
#include "utils/Logger.hpp"
#include <stb_image.h>
#include <imgui.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <algorithm>
#include <cmath>
#include <memory>
#include <mutex>

namespace optirad {

// ---------------------------------------------------------------------------
// Local shader-compilation helper
// ---------------------------------------------------------------------------
namespace {

static GLuint compileShader(GLenum type, const char* src) {
    GLuint sh = glCreateShader(type);
    glShaderSource(sh, 1, &src, nullptr);
    glCompileShader(sh);
    GLint ok = 0; glGetShaderiv(sh, GL_COMPILE_STATUS, &ok);
    if (!ok) {
        char buf[512]; glGetShaderInfoLog(sh, 512, nullptr, buf);
        Logger::error(std::string("View3D shader error: ") + buf);
    }
    return sh;
}

static GLuint linkProgram(const char* vsrc, const char* fsrc) {
    GLuint vs = compileShader(GL_VERTEX_SHADER, vsrc);
    GLuint fs = compileShader(GL_FRAGMENT_SHADER, fsrc);
    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs); glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok = 0; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char buf[512]; glGetProgramInfoLog(prog, 512, nullptr, buf);
        Logger::error(std::string("View3D program link error: ") + buf);
    }
    glDeleteShader(vs); glDeleteShader(fs);
    return prog;
}

static GLuint loadPNG(const char* path) {
    int w, h, ch;
    stbi_set_flip_vertically_on_load(true); // flip for OpenGL
    unsigned char* px = stbi_load(path, &w, &h, &ch, 4);
    if (!px) {
        Logger::error(std::string("View3D: failed to load image '") + path + "'");
        return 0;
    }
    GLuint tex = 0;
    glGenTextures(1, &tex);
    glBindTexture(GL_TEXTURE_2D, tex);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, px);
    glBindTexture(GL_TEXTURE_2D, 0);
    stbi_image_free(px);
    return tex;
}

} // anonymous namespace

struct View3D::Impl {
    std::unique_ptr<CubeRenderer> cubeRenderer;
    std::unique_ptr<AxisLabels> axisLabels;
    std::unique_ptr<VolumeRenderer> volumeRenderer;
    std::unique_ptr<StructureRenderer> structureRenderer;
    std::unique_ptr<BeamRenderer> beamRenderer;
    std::unique_ptr<PhaseSpaceRenderer> phaseSpaceRenderer;
    PatientData* lastLoadedData = nullptr;
};

View3D::View3D() : m_impl(std::make_unique<Impl>()) {}
View3D::~View3D() = default;

void View3D::init() {
    glEnable(GL_DEPTH_TEST);
    glLineWidth(2.0f);
    
    m_impl->cubeRenderer = std::make_unique<CubeRenderer>();
    m_impl->cubeRenderer->init();
    
    m_impl->axisLabels = std::make_unique<AxisLabels>();
    m_impl->axisLabels->init();
    
    m_impl->volumeRenderer = std::make_unique<VolumeRenderer>();
    m_impl->volumeRenderer->init();
    
    m_impl->structureRenderer = std::make_unique<StructureRenderer>();
    m_impl->structureRenderer->init();
    
    m_impl->beamRenderer = std::make_unique<BeamRenderer>();
    m_impl->beamRenderer->init();
    
    m_impl->phaseSpaceRenderer = std::make_unique<PhaseSpaceRenderer>();
    m_impl->phaseSpaceRenderer->init();
    
    initLogoRendering();
    initFBO();
}

void View3D::initFBO() {
    glGenFramebuffers(1, &m_fbo);
    glGenTextures(1, &m_fboTexture);
    glGenRenderbuffers(1, &m_fboDepthRBO);
}

void View3D::resizeFBO(int width, int height) {
    if (width <= 0 || height <= 0) return;
    if (width == m_fboWidth && height == m_fboHeight) return;
    
    m_fboWidth = width;
    m_fboHeight = height;
    
    // Color texture
    glBindTexture(GL_TEXTURE_2D, m_fboTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
    
    // Depth renderbuffer
    glBindRenderbuffer(GL_RENDERBUFFER, m_fboDepthRBO);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glBindRenderbuffer(GL_RENDERBUFFER, 0);
    
    // Attach to FBO
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_fboTexture, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_fboDepthRBO);
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void View3D::handleScroll(double yOffset) {
    std::lock_guard<std::mutex> lock(m_mouseMutex);
    m_scrollOffset = yOffset;
}

void View3D::handleMouseInput(GLFWwindow* window) {
    std::lock_guard<std::mutex> lock(m_mouseMutex);
    
    glfwGetFramebufferSize(window, &m_viewportWidth, &m_viewportHeight);
    
    double mouseX, mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);
    
    bool leftPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
    bool rightPressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS;
    bool middlePressed = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_MIDDLE) == GLFW_PRESS;
    
    if (ImGui::GetIO().WantCaptureMouse) {
        m_leftMousePressed = false;
        m_rightMousePressed = false;
        m_middleMousePressed = false;
        m_scrollOffset = 0.0;
        return;
    }
    
    double deltaX = mouseX - m_lastMouseX;
    double deltaY = mouseY - m_lastMouseY;
    
    if (leftPressed && m_leftMousePressed) {
        const float sensitivity = 0.005f;
        m_yaw -= static_cast<float>(deltaX) * sensitivity;
        m_pitch += static_cast<float>(deltaY) * sensitivity;
        const float maxPitch = glm::half_pi<float>() - 0.01f;
        m_pitch = std::clamp(m_pitch, -maxPitch, maxPitch);
    }
    
    if (rightPressed && m_rightMousePressed) {
        const float panSpeed = 0.002f * m_cameraDistance;
        glm::vec3 forward(std::sin(m_yaw), 0.0f, std::cos(m_yaw));
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        m_target -= right * static_cast<float>(deltaX) * panSpeed;
        m_target += up * static_cast<float>(deltaY) * panSpeed;
    }
    
    if (m_scrollOffset != 0.0) {
        const float zoomSpeed = 0.1f;
        m_cameraDistance -= static_cast<float>(m_scrollOffset) * zoomSpeed * m_cameraDistance;
        m_cameraDistance = std::clamp(m_cameraDistance, 0.01f, 50.0f);
        m_scrollOffset = 0.0;
    }
    
    m_lastMouseX = mouseX;
    m_lastMouseY = mouseY;
    m_leftMousePressed = leftPressed;
    m_rightMousePressed = rightPressed;
    m_middleMousePressed = middlePressed;
}

void View3D::setPatientData(PatientData* data) {
    if (!m_impl->volumeRenderer) return;
    m_impl->volumeRenderer->setPatientData(data);
    m_impl->structureRenderer->setPatientData(data);
    
    // Adjust camera to fit volume when data is first loaded (only once)
    if (data && data->getCTVolume() && data != m_impl->lastLoadedData) {
        m_cameraDistance = 3.0f;
        m_impl->lastLoadedData = data;
    }
}

void View3D::setStf(const Stf* stf) {
    if (!m_impl->beamRenderer) return;
    m_impl->beamRenderer->setStf(stf);
}

BeamRenderer* View3D::getBeamRenderer() {
    return m_impl->beamRenderer.get();
}

void View3D::setPhaseSpaceSources(const std::vector<const PhaseSpaceBeamSource*>& sources) {
    if (!m_impl->phaseSpaceRenderer) return;
    m_impl->phaseSpaceRenderer->setSources(sources);
}

PhaseSpaceRenderer* View3D::getPhaseSpaceRenderer() {
    return m_impl->phaseSpaceRenderer.get();
}

void View3D::handleImGuiInput() {
    // This is called when the 3D View ImGui window is hovered.
    // We read ImGui IO state instead of raw GLFW.
    ImGuiIO& io = ImGui::GetIO();
    
    ImVec2 mousePos = io.MousePos;
    double mouseX = static_cast<double>(mousePos.x);
    double mouseY = static_cast<double>(mousePos.y);
    
    bool leftPressed = ImGui::IsMouseDown(ImGuiMouseButton_Left);
    bool rightPressed = ImGui::IsMouseDown(ImGuiMouseButton_Right);
    
    double deltaX = mouseX - m_lastMouseX;
    double deltaY = mouseY - m_lastMouseY;
    
    if (leftPressed && m_leftMousePressed) {
        const float sensitivity = 0.005f;
        m_yaw -= static_cast<float>(deltaX) * sensitivity;
        m_pitch += static_cast<float>(deltaY) * sensitivity;
        const float maxPitch = glm::half_pi<float>() - 0.01f;
        m_pitch = std::clamp(m_pitch, -maxPitch, maxPitch);
    }
    
    if (rightPressed && m_rightMousePressed) {
        const float panSpeed = 0.002f * m_cameraDistance;
        glm::vec3 forward(std::sin(m_yaw), 0.0f, std::cos(m_yaw));
        glm::vec3 right = glm::normalize(glm::cross(forward, glm::vec3(0.0f, 1.0f, 0.0f)));
        glm::vec3 up(0.0f, 1.0f, 0.0f);
        m_target -= right * static_cast<float>(deltaX) * panSpeed;
        m_target += up * static_cast<float>(deltaY) * panSpeed;
    }
    
    // Scroll zoom
    float scroll = io.MouseWheel;
    if (scroll != 0.0f) {
        const float zoomSpeed = 0.1f;
        m_cameraDistance -= scroll * zoomSpeed * m_cameraDistance;
        m_cameraDistance = std::clamp(m_cameraDistance, 0.01f, 50.0f);
    }
    
    m_lastMouseX = mouseX;
    m_lastMouseY = mouseY;
    m_leftMousePressed = leftPressed;
    m_rightMousePressed = rightPressed;
}

void View3D::render(int width, int height) {
    if (width <= 0 || height <= 0) return;
    
    m_viewportWidth = width;
    m_viewportHeight = height;
    
    // Ensure FBO is correct size
    resizeFBO(width, height);
    
    // Bind FBO for off-screen rendering
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
    glViewport(0, 0, width, height);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if (m_impl->lastLoadedData == nullptr) {
        // ---- No patient loaded: show logo splash ----
        renderLogo(width, height);
    } else {
        // ---- Patient loaded: render 3D scene ----
        if (m_darkMode) {
            // Dark mode: solid black background
            glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        } else {
            // Light mode: white-to-sky gradient background, then depth-only clear
            glClear(GL_DEPTH_BUFFER_BIT);
            glDisable(GL_DEPTH_TEST);
            if (m_gradShader && m_gradVAO) {
                glUseProgram(m_gradShader);
                glBindVertexArray(m_gradVAO);
                glDrawArrays(GL_TRIANGLES, 0, 6);
                glBindVertexArray(0);
                glUseProgram(0);
            }
        }
        glEnable(GL_DEPTH_TEST);

        float cosPitch = std::cos(m_pitch), sinPitch = std::sin(m_pitch);
        float cosYaw = std::cos(m_yaw), sinYaw = std::sin(m_yaw);
        
        glm::vec3 cameraOffset(
            m_cameraDistance * cosPitch * sinYaw,
            m_cameraDistance * sinPitch,
            m_cameraDistance * cosPitch * cosYaw
        );
        
        glm::vec3 cameraPos = m_target + cameraOffset;
        glm::mat4 view = glm::lookAt(cameraPos, m_target, glm::vec3(0.0f, 1.0f, 0.0f));
        
        float aspect = static_cast<float>(width) / static_cast<float>(std::max(height, 1));
        glm::mat4 projection = glm::perspective(glm::radians(45.0f), aspect, 0.01f, 100.0f);
        
        // Render structures (semi-transparent)
        m_impl->structureRenderer->render(view, projection);
        
        // Render beams / rays / isocenter
        m_impl->beamRenderer->render(view, projection, cameraPos);
        
        // Render phase-space particles
        m_impl->phaseSpaceRenderer->render(view, projection, cameraPos);

        // Render orientation cube and labels on top
        m_impl->cubeRenderer->render(view, projection);
        m_impl->axisLabels->render(view, projection);
    }
    
    // Unbind FBO - go back to default framebuffer
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

// ---------------------------------------------------------------------------
// Logo rendering (no-patient splash screen)
// ---------------------------------------------------------------------------

void View3D::initLogoRendering() {
    // ---- Logo shader (textured quad, uniform scale) ----
    const char* logoVS = R"(
        #version 330 core
        layout(location = 0) in vec2 aPos;
        layout(location = 1) in vec2 aUV;
        uniform vec2 uScale;
        out vec2 vUV;
        void main() {
            gl_Position = vec4(aPos.x * uScale.x, aPos.y * uScale.y, 0.0, 1.0);
            vUV = aUV;
        }
    )";
    const char* logoFS = R"(
        #version 330 core
        in vec2 vUV;
        uniform sampler2D uTex;
        out vec4 FragColor;
        void main() {
            FragColor = texture(uTex, vUV);
        }
    )";
    m_logoShader = linkProgram(logoVS, logoFS);

    // Unit quad (local space [-1,1] x [-1,1]) — two triangles
    float logoVerts[] = {
        // x      y     u     v
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f, -1.0f,  1.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
         1.0f,  1.0f,  1.0f, 1.0f,
        -1.0f,  1.0f,  0.0f, 1.0f,
    };
    glGenVertexArrays(1, &m_logoVAO);
    glGenBuffers(1, &m_logoVBO);
    glBindVertexArray(m_logoVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_logoVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(logoVerts), logoVerts, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);

    // ---- Gradient shader (for light mode background) ----
    const char* gradVS = R"(
        #version 330 core
        layout(location = 0) in vec2 aPos;
        layout(location = 1) in vec4 aColor;
        out vec4 vColor;
        void main() {
            gl_Position = vec4(aPos, 0.0, 1.0);
            vColor = aColor;
        }
    )";
    const char* gradFS = R"(
        #version 330 core
        in vec4 vColor;
        out vec4 FragColor;
        void main() {
            FragColor = vColor;
        }
    )";
    m_gradShader = linkProgram(gradVS, gradFS);

    // Sky blue #5EC2EC = (94/255, 194/255, 236/255)
    constexpr float SR = 94.0f / 255.0f;
    constexpr float SG = 194.0f / 255.0f;
    constexpr float SB = 236.0f / 255.0f;
    // Full-viewport quad. Top = white (y=+1), Bottom = sky (y=-1)
    float gradVerts[] = {
        // x      y     r     g     b     a
        -1.0f, -1.0f,   SR,   SG,   SB,  1.0f, // bottom-left  sky
         1.0f, -1.0f,   SR,   SG,   SB,  1.0f, // bottom-right sky
         1.0f,  1.0f,  1.0f, 1.0f, 1.0f, 1.0f, // top-right    white
        -1.0f, -1.0f,   SR,   SG,   SB,  1.0f, // bottom-left  sky
         1.0f,  1.0f,  1.0f, 1.0f, 1.0f, 1.0f, // top-right    white
        -1.0f,  1.0f,  1.0f, 1.0f, 1.0f, 1.0f, // top-left     white
    };
    glGenVertexArrays(1, &m_gradVAO);
    glGenBuffers(1, &m_gradVBO);
    glBindVertexArray(m_gradVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_gradVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(gradVerts), gradVerts, GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(2 * sizeof(float)));
    glBindVertexArray(0);

    // ---- Load logo textures ----
    m_logoTexDark  = loadPNG(OPTIRAD_GUI_IMG_DIR "/LaTIM_white_black.png");
    m_logoTexLight = loadPNG(OPTIRAD_GUI_IMG_DIR "/LaTIM.png");
}

void View3D::renderLogo(int width, int height) {
    // Determine logo pixel size: 40% of the shorter viewport dimension
    float logoPx = 0.40f * static_cast<float>(std::min(width, height));
    // Convert to NDC half-extents (NDC range is [-1,1] = 2 units per axis)
    float scaleX = logoPx / static_cast<float>(width);
    float scaleY = logoPx / static_cast<float>(height);

    glDisable(GL_DEPTH_TEST);

    if (m_darkMode) {
        // Black background
        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);
    } else {
        // White → sky gradient background
        glUseProgram(m_gradShader);
        glBindVertexArray(m_gradVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
    }

    // Render logo on top with alpha blending
    GLuint logoTex = m_darkMode ? m_logoTexDark : m_logoTexLight;
    if (logoTex && m_logoShader && m_logoVAO) {
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

        glUseProgram(m_logoShader);
        glUniform2f(glGetUniformLocation(m_logoShader, "uScale"), scaleX, scaleY);
        glUniform1i(glGetUniformLocation(m_logoShader, "uTex"), 0);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, logoTex);
        glBindVertexArray(m_logoVAO);
        glDrawArrays(GL_TRIANGLES, 0, 6);
        glBindVertexArray(0);
        glBindTexture(GL_TEXTURE_2D, 0);

        glDisable(GL_BLEND);
    }

    glEnable(GL_DEPTH_TEST);
    glUseProgram(0);
}

void View3D::cleanupLogoRendering() {
    if (m_logoShader)   { glDeleteProgram(m_logoShader);  m_logoShader  = 0; }
    if (m_gradShader)   { glDeleteProgram(m_gradShader);  m_gradShader  = 0; }
    if (m_logoVAO)      { glDeleteVertexArrays(1, &m_logoVAO); m_logoVAO = 0; }
    if (m_logoVBO)      { glDeleteBuffers(1, &m_logoVBO);      m_logoVBO = 0; }
    if (m_gradVAO)      { glDeleteVertexArrays(1, &m_gradVAO); m_gradVAO = 0; }
    if (m_gradVBO)      { glDeleteBuffers(1, &m_gradVBO);      m_gradVBO = 0; }
    if (m_logoTexDark)  { glDeleteTextures(1, &m_logoTexDark);  m_logoTexDark  = 0; }
    if (m_logoTexLight) { glDeleteTextures(1, &m_logoTexLight); m_logoTexLight = 0; }
}

// ---------------------------------------------------------------------------

void View3D::cleanup() {
    cleanupLogoRendering();
    if (m_fbo) { glDeleteFramebuffers(1, &m_fbo); m_fbo = 0; }
    if (m_fboTexture) { glDeleteTextures(1, &m_fboTexture); m_fboTexture = 0; }
    if (m_fboDepthRBO) { glDeleteRenderbuffers(1, &m_fboDepthRBO); m_fboDepthRBO = 0; }
    m_impl->cubeRenderer->cleanup();
    m_impl->axisLabels->cleanup();
    m_impl->volumeRenderer->cleanup();
    m_impl->structureRenderer->cleanup();
    m_impl->beamRenderer->cleanup();
    m_impl->phaseSpaceRenderer->cleanup();
}
}