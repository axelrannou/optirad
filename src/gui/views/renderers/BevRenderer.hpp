#pragma once

#include "StructureRenderer.hpp"
#include "core/Aperture.hpp"
#include "core/Beam.hpp"
#include <vector>

namespace optirad {

/// Source-perspective anatomy and MLC plane, in the same world frame as View3D.
class BevRenderer {
public:
    ~BevRenderer();

    GLuint render(PatientData* patient, const Beam& beam,
                  const std::vector<const Aperture*>& apertures,
                  const std::vector<double>& boundaries,
                  double extentX, double extentZ, int width, int height,
                  float zoom);

private:
    void init();
    void resize(int width, int height);

    StructureRenderer m_structures;
    GLuint m_fbo = 0, m_texture = 0, m_depth = 0;
    GLuint m_program = 0, m_vao = 0, m_vbo = 0;
    int m_width = 0, m_height = 0;
};

} // namespace optirad
