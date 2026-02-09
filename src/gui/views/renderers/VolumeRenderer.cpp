#include "VolumeRenderer.hpp"
#include "core/PatientData.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <vector>
#include <algorithm>
#include <cmath>

namespace optirad {

VolumeRenderer::VolumeRenderer() = default;
VolumeRenderer::~VolumeRenderer() = default;

void VolumeRenderer::init() {
    createProxyCube();
    createShaders();

    // Full-screen quad for final raycast pass
    float quadVerts[] = {
        -1.f, -1.f, 0.f, 0.f,
         1.f, -1.f, 1.f, 0.f,
         1.f,  1.f, 1.f, 1.f,
        -1.f, -1.f, 0.f, 0.f,
         1.f,  1.f, 1.f, 1.f,
        -1.f,  1.f, 0.f, 1.f,
    };
    glGenVertexArrays(1, &m_quadVAO);
    glGenBuffers(1, &m_quadVBO);
    glBindVertexArray(m_quadVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_quadVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(2 * sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);
}

void VolumeRenderer::createProxyCube() {
    // Unit cube [0,1]^3 used as proxy geometry
    float verts[] = {
        0,0,0,  1,0,0,  1,1,0,  0,1,0,
        0,0,1,  1,0,1,  1,1,1,  0,1,1,
    };
    unsigned int indices[] = {
        // front
        0,1,2, 2,3,0,
        // right
        1,5,6, 6,2,1,
        // back
        5,4,7, 7,6,5,
        // left
        4,0,3, 3,7,4,
        // top
        3,2,6, 6,7,3,
        // bottom
        4,5,1, 1,0,4,
    };

    glGenVertexArrays(1, &m_cubeVAO);
    glGenBuffers(1, &m_cubeVBO);
    glGenBuffers(1, &m_cubeEBO);

    glBindVertexArray(m_cubeVAO);
    glBindBuffer(GL_ARRAY_BUFFER, m_cubeVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, m_cubeEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

void VolumeRenderer::createShaders() {
    // --- Position shader: renders cube vertex positions as colors ---
    const char* posVS = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        out vec3 vPos;
        uniform mat4 MVP;
        void main() {
            vPos = aPos;
            gl_Position = MVP * vec4(aPos, 1.0);
        }
    )";
    const char* posFS = R"(
        #version 330 core
        in vec3 vPos;
        out vec4 FragColor;
        void main() {
            FragColor = vec4(vPos, 1.0);
        }
    )";

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &posVS, nullptr);
    glCompileShader(vs);
    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &posFS, nullptr);
    glCompileShader(fs);
    m_positionShader = glCreateProgram();
    glAttachShader(m_positionShader, vs);
    glAttachShader(m_positionShader, fs);
    glLinkProgram(m_positionShader);
    glDeleteShader(vs);
    glDeleteShader(fs);

    // --- Raycast shader: uses front and back face textures ---
    const char* rayVS = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        out vec3 vPos;
        uniform mat4 MVP;
        void main() {
            vPos = aPos;
            gl_Position = MVP * vec4(aPos, 1.0);
        }
    )";
    const char* rayFS = R"(
        #version 330 core
        in vec3 vPos;
        out vec4 FragColor;

        uniform sampler3D volumeTex;
        uniform sampler2D backFaceTex;
        uniform float windowMin;
        uniform float windowMax;
        uniform int numSteps;

        void main() {
            // Entry point is the front face (vPos is [0,1])
            vec3 entryPoint = vPos;
            
            // Exit point from back face texture (need to convert screen coords)
            vec2 screenCoord = gl_FragCoord.xy;
            ivec2 vpSize = textureSize(backFaceTex, 0);
            vec2 texCoord = screenCoord / vec2(vpSize);
            vec3 exitPoint = texture(backFaceTex, texCoord).xyz;
            
            // If no valid back face, skip
            if (exitPoint == vec3(0.0)) {
                discard;
            }

            vec3 rayDir = exitPoint - entryPoint;
            float rayLen = length(rayDir);
            
            if (rayLen < 0.001) {
                discard;
            }
            
            rayDir = normalize(rayDir);
            float stepSize = rayLen / float(numSteps);

            vec4 accum = vec4(0.0);
            vec3 pos = entryPoint;

            for (int i = 0; i < numSteps; ++i) {
                // Skip if outside [0,1]
                if (any(lessThan(pos, vec3(0.0))) || any(greaterThan(pos, vec3(1.0)))) {
                    pos += rayDir * stepSize;
                    continue;
                }

                float sample = texture(volumeTex, pos).r;

                // Map to [0,1] using window/level
                float intensity = clamp((sample - windowMin) / (windowMax - windowMin), 0.0, 1.0);

                // Transfer function
                float alpha = 0.0;
                if (intensity > 0.15) {
                    alpha = intensity * 0.3;
                }

                vec3 color = vec3(intensity);

                // Front-to-back compositing
                accum.rgb += (1.0 - accum.a) * alpha * color;
                accum.a   += (1.0 - accum.a) * alpha;

                if (accum.a > 0.95) break;
                
                pos += rayDir * stepSize;
            }

            FragColor = accum;
        }
    )";

    vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &rayVS, nullptr);
    glCompileShader(vs);
    fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &rayFS, nullptr);
    glCompileShader(fs);
    m_raycastShader = glCreateProgram();
    glAttachShader(m_raycastShader, vs);
    glAttachShader(m_raycastShader, fs);
    glLinkProgram(m_raycastShader);
    glDeleteShader(vs);
    glDeleteShader(fs);
}

void VolumeRenderer::createBackFBO(int width, int height) {
    if (m_backFBO) {
        glDeleteFramebuffers(1, &m_backFBO);
        glDeleteTextures(1, &m_backTexture);
        glDeleteRenderbuffers(1, &m_backDepthRBO);
    }

    glGenFramebuffers(1, &m_backFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_backFBO);

    glGenTextures(1, &m_backTexture);
    glBindTexture(GL_TEXTURE_2D, m_backTexture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_backTexture, 0);

    glGenRenderbuffers(1, &m_backDepthRBO);
    glBindRenderbuffer(GL_RENDERBUFFER, m_backDepthRBO);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_backDepthRBO);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    m_lastWidth = width;
    m_lastHeight = height;
}

void VolumeRenderer::setPatientData(PatientData* data) {
    if (m_patientData == data) return;
    m_patientData = data;
    m_needsUpload = true;
}

void VolumeRenderer::setWindowLevel(int windowWidth, int windowCenter) {
    m_windowWidth = windowWidth;
    m_windowCenter = windowCenter;
}

void VolumeRenderer::uploadVolumeTexture() {
    if (!m_patientData || !m_patientData->getCTVolume()) return;

    auto* ct = m_patientData->getCTVolume();
    auto dims = ct->getGrid().getDimensions();
    auto spacing = ct->getGrid().getSpacing();

    m_dimX = static_cast<int>(dims[0]);
    m_dimY = static_cast<int>(dims[1]);
    m_dimZ = static_cast<int>(dims[2]);

    // Compute physical size and normalize so the largest axis = 1.0
    float physX = m_dimX * static_cast<float>(spacing[0]);
    float physY = m_dimY * static_cast<float>(spacing[1]);
    float physZ = m_dimZ * static_cast<float>(spacing[2]);
    float maxPhys = std::max({physX, physY, physZ});
    m_volumeScale = glm::vec3(physX / maxPhys, physY / maxPhys, physZ / maxPhys);

    // Normalize HU values to [0,1] range for 16-bit texture
    // HU range: roughly -1024 to 3071
    size_t totalVoxels = static_cast<size_t>(m_dimX) * m_dimY * m_dimZ;
    std::vector<uint16_t> normalized(totalVoxels);

    const int16_t* src = ct->data();
    for (size_t i = 0; i < totalVoxels; ++i) {
        float val = (static_cast<float>(src[i]) + 1024.0f) / 4096.0f;
        val = std::clamp(val, 0.0f, 1.0f);
        normalized[i] = static_cast<uint16_t>(val * 65535.0f);
    }

    if (!m_volumeTexture) glGenTextures(1, &m_volumeTexture);
    glBindTexture(GL_TEXTURE_3D, m_volumeTexture);
    glTexImage3D(GL_TEXTURE_3D, 0, GL_R16, m_dimX, m_dimY, m_dimZ, 0,
                 GL_RED, GL_UNSIGNED_SHORT, normalized.data());
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    m_needsUpload = false;
}

void VolumeRenderer::render(const glm::mat4& view, const glm::mat4& projection, const glm::vec3& cameraPos) {
    if (!m_patientData || !m_patientData->getCTVolume()) return;

    if (m_needsUpload) uploadVolumeTexture();
    if (!m_volumeTexture) return;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int vpW = viewport[2], vpH = viewport[3];

    if (vpW != m_lastWidth || vpH != m_lastHeight) {
        createBackFBO(vpW, vpH);
    }

    // Model matrix: scale and center
    glm::mat4 model = glm::mat4(1.0f);
    model = glm::scale(model, m_volumeScale);
    model = glm::translate(model, glm::vec3(-0.5f));

    glm::mat4 MVP = projection * view * model;

    // --- Pass 1: Render back faces to texture ---
    glBindFramebuffer(GL_FRAMEBUFFER, m_backFBO);
    glClearColor(0, 0, 0, 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_CULL_FACE);
    glCullFace(GL_FRONT); // Cull front, draw back
    glEnable(GL_DEPTH_TEST);

    glUseProgram(m_positionShader);
    glUniformMatrix4fv(glGetUniformLocation(m_positionShader, "MVP"), 1, GL_FALSE, glm::value_ptr(MVP));
    glBindVertexArray(m_cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

    // --- Pass 2: Render front faces with raycasting ---
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glCullFace(GL_BACK); // Cull back, draw front
    glDisable(GL_DEPTH_TEST); // Volume renders on top

    float windowMin = (static_cast<float>(m_windowCenter - m_windowWidth / 2) + 1024.0f) / 4096.0f;
    float windowMax = (static_cast<float>(m_windowCenter + m_windowWidth / 2) + 1024.0f) / 4096.0f;

    glUseProgram(m_raycastShader);
    glUniform1f(glGetUniformLocation(m_raycastShader, "windowMin"), windowMin);
    glUniform1f(glGetUniformLocation(m_raycastShader, "windowMax"), windowMax);
    glUniform1i(glGetUniformLocation(m_raycastShader, "numSteps"), 512);
    glUniformMatrix4fv(glGetUniformLocation(m_raycastShader, "MVP"), 1, GL_FALSE, glm::value_ptr(MVP));

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_3D, m_volumeTexture);
    glUniform1i(glGetUniformLocation(m_raycastShader, "volumeTex"), 0);

    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, m_backTexture);
    glUniform1i(glGetUniformLocation(m_raycastShader, "backFaceTex"), 1);

    // Draw front faces (each fragment gets entry point as vPos)
    glBindVertexArray(m_cubeVAO);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);

    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glActiveTexture(GL_TEXTURE0);
    glUseProgram(0);
}

void VolumeRenderer::cleanup() {
    if (m_cubeVAO) { glDeleteVertexArrays(1, &m_cubeVAO); m_cubeVAO = 0; }
    if (m_cubeVBO) { glDeleteBuffers(1, &m_cubeVBO); m_cubeVBO = 0; }
    if (m_cubeEBO) { glDeleteBuffers(1, &m_cubeEBO); m_cubeEBO = 0; }
    if (m_quadVAO) { glDeleteVertexArrays(1, &m_quadVAO); m_quadVAO = 0; }
    if (m_quadVBO) { glDeleteBuffers(1, &m_quadVBO); m_quadVBO = 0; }
    if (m_volumeTexture) { glDeleteTextures(1, &m_volumeTexture); m_volumeTexture = 0; }
    if (m_backFBO) { glDeleteFramebuffers(1, &m_backFBO); m_backFBO = 0; }
    if (m_backTexture) { glDeleteTextures(1, &m_backTexture); m_backTexture = 0; }
    if (m_backDepthRBO) { glDeleteRenderbuffers(1, &m_backDepthRBO); m_backDepthRBO = 0; }
    if (m_positionShader) { glDeleteProgram(m_positionShader); m_positionShader = 0; }
    if (m_raycastShader) { glDeleteProgram(m_raycastShader); m_raycastShader = 0; }
}

}
