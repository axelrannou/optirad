#include "BevRenderer.hpp"
#include "utils/Logger.hpp"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

namespace optirad {
namespace {

GLuint compile(GLenum type, const char* source) {
    GLuint shader = glCreateShader(type);
    glShaderSource(shader, 1, &source, nullptr);
    glCompileShader(shader);
    GLint success = 0;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
    if (!success) {
        char log[1024];
        glGetShaderInfoLog(shader, sizeof(log), nullptr, log);
        glDeleteShader(shader);
        Logger::error(std::string("BEV shader: ") + log);
        throw std::runtime_error("Failed to compile BEV shader");
    }
    return shader;
}

} // namespace

BevRenderer::~BevRenderer() {
    m_structures.cleanup();
    if (m_vbo) glDeleteBuffers(1, &m_vbo);
    if (m_vao) glDeleteVertexArrays(1, &m_vao);
    if (m_program) glDeleteProgram(m_program);
    if (m_depth) glDeleteRenderbuffers(1, &m_depth);
    if (m_texture) glDeleteTextures(1, &m_texture);
    if (m_fbo) glDeleteFramebuffers(1, &m_fbo);
}

void BevRenderer::init() {
    const char* vertex = R"(
        #version 330 core
        layout(location = 0) in vec3 position;
        uniform mat4 mvp;
        void main() { gl_Position = mvp * vec4(position, 1.0); }
    )";
    const char* fragment = R"(
        #version 330 core
        uniform vec4 color;
        out vec4 result;
        void main() { result = color; }
    )";
    GLuint vs = compile(GL_VERTEX_SHADER, vertex);
    GLuint fs = compile(GL_FRAGMENT_SHADER, fragment);
    m_program = glCreateProgram();
    glAttachShader(m_program, vs);
    glAttachShader(m_program, fs);
    glLinkProgram(m_program);
    glDeleteShader(vs);
    glDeleteShader(fs);
    GLint success = 0;
    glGetProgramiv(m_program, GL_LINK_STATUS, &success);
    if (!success) {
        char log[1024];
        glGetProgramInfoLog(m_program, sizeof(log), nullptr, log);
        Logger::error(std::string("BEV program: ") + log);
        throw std::runtime_error("Failed to link BEV shader");
    }
    glGenFramebuffers(1, &m_fbo);
    glGenTextures(1, &m_texture);
    glGenRenderbuffers(1, &m_depth);
    glGenVertexArrays(1, &m_vao);
    glGenBuffers(1, &m_vbo);
    glBindVertexArray(m_vao);
    glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(glm::vec3), nullptr);
    glBindVertexArray(0);
    m_structures.init();
}

void BevRenderer::resize(int width, int height) {
    glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
    if (width == m_width && height == m_height) return;
    glBindTexture(GL_TEXTURE_2D, m_texture);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindRenderbuffer(GL_RENDERBUFFER, m_depth);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_texture, 0);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_depth);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE) {
        Logger::error("BEV framebuffer is incomplete");
        throw std::runtime_error("Failed to create BEV framebuffer");
    }
    m_width = width;
    m_height = height;
}

GLuint BevRenderer::render(PatientData* patient, const Beam& beam,
        const std::vector<const Aperture*>& apertures, const std::vector<double>& boundaries,
        double extentX, double extentZ, int width, int height,
        float zoom) {
    if (width <= 0 || height <= 0 || !std::isfinite(beam.getSAD()) ||
        !std::isfinite(beam.getSCD()) || beam.getSAD() <= 0.0 ||
        beam.getSCD() <= 0.0 || beam.getSCD() >= beam.getSAD()) {
        Logger::error("BEV requires a positive viewport and 0 < SCD < SAD");
        return 0;
    }
    GLint drawFbo, readFbo, viewport[4], program, vao, arrayBuffer, renderbuffer, activeTexture;
    GLint textures[2], blendSrcRgb, blendDstRgb, blendSrcAlpha, blendDstAlpha, blendEqRgb, blendEqAlpha;
    GLfloat clearColor[4];
    glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &drawFbo);
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &readFbo);
    glGetIntegerv(GL_VIEWPORT, viewport);
    glGetIntegerv(GL_CURRENT_PROGRAM, &program);
    glGetIntegerv(GL_VERTEX_ARRAY_BINDING, &vao);
    glGetIntegerv(GL_ARRAY_BUFFER_BINDING, &arrayBuffer);
    glGetIntegerv(GL_RENDERBUFFER_BINDING, &renderbuffer);
    glGetIntegerv(GL_ACTIVE_TEXTURE, &activeTexture);
    for (int unit = 0; unit < 2; ++unit) {
        glActiveTexture(GL_TEXTURE0 + unit);
        glGetIntegerv(GL_TEXTURE_BINDING_2D, &textures[unit]);
    }
    glActiveTexture(GL_TEXTURE0);
    glGetIntegerv(GL_BLEND_SRC_RGB, &blendSrcRgb);
    glGetIntegerv(GL_BLEND_DST_RGB, &blendDstRgb);
    glGetIntegerv(GL_BLEND_SRC_ALPHA, &blendSrcAlpha);
    glGetIntegerv(GL_BLEND_DST_ALPHA, &blendDstAlpha);
    glGetIntegerv(GL_BLEND_EQUATION_RGB, &blendEqRgb);
    glGetIntegerv(GL_BLEND_EQUATION_ALPHA, &blendEqAlpha);
    glGetFloatv(GL_COLOR_CLEAR_VALUE, clearColor);
    const GLboolean scissor = glIsEnabled(GL_SCISSOR_TEST);
    const GLboolean depth = glIsEnabled(GL_DEPTH_TEST);
    const GLboolean cull = glIsEnabled(GL_CULL_FACE);
    const GLboolean blend = glIsEnabled(GL_BLEND);
    GLboolean depthMask;
    glGetBooleanv(GL_DEPTH_WRITEMASK, &depthMask);

    if (!m_fbo) init();
    resize(width, height);
    glViewport(0, 0, width, height);
    glDisable(GL_SCISSOR_TEST);
    glDisable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE);
    glClearColor(0.055f, 0.07f, 0.085f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // StructureRenderer uses metres and a -90 degree rotation around LPS X.
    const auto rotation = getRotationMatrix(beam.getGantryAngle(), beam.getCouchAngle());
    auto world = [&](const Vec3& local) {
        const auto lps = vecAdd(beam.getIsocenter(), rotation * local);
        return glm::vec3(lps[0], lps[2], -lps[1]) * 0.001f;
    };
    const float sad = static_cast<float>(beam.getSAD());
    const float aspect = static_cast<float>(width) / height;
    const float halfHeight = static_cast<float>(std::max(extentZ, extentX / aspect)) * 1.15f / zoom;
    const float nearPlane = 0.001f;
    const float halfNear = halfHeight / sad * nearPlane;
    const auto target = world({0.0, 0.0, 0.0});
    const auto camera = world({0.0, -sad, 0.0});
    const auto up = glm::normalize(world({0.0, 0.0, 1.0}) - world({0.0, 0.0, 0.0}));
    const auto view = glm::lookAt(camera, target, up);
    const auto projection = glm::frustum(-halfNear * aspect, halfNear * aspect,
                                         -halfNear, halfNear, nearPlane, sad * 0.001f + 10.0f);
    m_structures.setPatientData(patient);
    m_structures.render(view, projection);

    // Leaf positions/widths are specified at isocenter, not at the MLC plane.
    const double magnification = beam.getSCD() / beam.getSAD();
    auto leafPoint = [&](double x, double z) {
        return world({x * magnification, beam.getSCD() - beam.getSAD(), z * magnification});
    };
    std::vector<glm::vec3> faces, edges, tips;
    auto rectangle = [&](double left, double right, double bottom, double top, bool filled) {
        const glm::vec3 a = leafPoint(left, bottom), b = leafPoint(right, bottom);
        const glm::vec3 c = leafPoint(right, top), d = leafPoint(left, top);
        if (filled) faces.insert(faces.end(), {a, b, c, a, c, d});
        edges.insert(edges.end(), {a, b, b, c, c, d, d, a});
    };
    for (const auto* aperture : apertures) {
        if (!aperture || aperture->bankA.size() != aperture->bankB.size()) {
            Logger::error("BEV aperture requires matching leaf banks");
            continue;
        }
        const size_t pairs = aperture->bankA.size();
        if (boundaries.size() != pairs + 1) {
            Logger::error("BEV leaf boundaries do not match aperture banks");
            continue;
        }
        for (size_t pair = 0; pair < pairs; ++pair) {
            const double a = aperture->bankA[pair], b = aperture->bankB[pair];
            const double bottom = boundaries[pair], top = boundaries[pair + 1];
            if (apertures.size() == 1) {
                rectangle(std::min(-extentX, a), a, bottom, top, true);
                rectangle(b, std::max(extentX, b), bottom, top, true);
            } else {
                rectangle(a, b, bottom, top, false);
            }
            tips.insert(tips.end(), {leafPoint(a, bottom), leafPoint(a, top),
                                    leafPoint(b, bottom), leafPoint(b, top)});
        }
    }
    glDisable(GL_CULL_FACE);
    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendEquation(GL_FUNC_ADD);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glUseProgram(m_program);
    const auto mvp = projection * view;
    glUniformMatrix4fv(glGetUniformLocation(m_program, "mvp"), 1, GL_FALSE, glm::value_ptr(mvp));
    glBindVertexArray(m_vao);
    auto draw = [&](const std::vector<glm::vec3>& vertices, GLenum mode, const glm::vec4& color) {
        glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(glm::vec3), vertices.data(), GL_STREAM_DRAW);
        glUniform4fv(glGetUniformLocation(m_program, "color"), 1, glm::value_ptr(color));
        glDrawArrays(mode, 0, static_cast<GLsizei>(vertices.size()));
    };
    draw(faces, GL_TRIANGLES, {0.35f, 0.42f, 0.49f, 0.65f});
    draw(edges, GL_LINES, {0.65f, 0.72f, 0.8f, 0.7f});
    draw(tips, GL_LINES, {0.2f, 0.95f, 0.75f, 1.0f});
    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, arrayBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, renderbuffer);
    glUseProgram(program);
    for (int unit = 0; unit < 2; ++unit) {
        glActiveTexture(GL_TEXTURE0 + unit);
        glBindTexture(GL_TEXTURE_2D, textures[unit]);
    }
    glActiveTexture(activeTexture);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, drawFbo);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, readFbo);
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    if (scissor) glEnable(GL_SCISSOR_TEST);
    if (depth) glEnable(GL_DEPTH_TEST); else glDisable(GL_DEPTH_TEST);
    if (cull) glEnable(GL_CULL_FACE); else glDisable(GL_CULL_FACE);
    if (!blend) glDisable(GL_BLEND);
    glBlendFuncSeparate(blendSrcRgb, blendDstRgb, blendSrcAlpha, blendDstAlpha);
    glBlendEquationSeparate(blendEqRgb, blendEqAlpha);
    glDepthMask(depthMask);
    glClearColor(clearColor[0], clearColor[1], clearColor[2], clearColor[3]);
    return m_texture;
}

} // namespace optirad
