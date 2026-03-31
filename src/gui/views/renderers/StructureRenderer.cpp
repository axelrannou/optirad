#include "StructureRenderer.hpp"
#include "MarchingCubesTables.hpp"
#include "core/PatientData.hpp"
#include <glm/gtc/type_ptr.hpp>
#include <glm/glm.hpp>
#include <algorithm>
#include <cmath>
#include <unordered_map>
#include <functional>
#include "utils/Logger.hpp"

namespace optirad {

StructureRenderer::StructureRenderer() = default;
StructureRenderer::~StructureRenderer() = default;

void StructureRenderer::init() {
    // ── Geometry pass shader (writes weighted accumulation buffers) ──────────
    // WBOIT accumulation pass: outputs (weightedColor, weightedAlpha) into
    // two MRT targets. No depth writes — all fragments contribute.
    const char* accumVertSrc = R"(
        #version 330 core
        layout (location = 0) in vec3 aPos;
        layout (location = 1) in vec3 aNormal;

        out vec3 FragPos;
        out vec3 Normal;
        out vec3 ViewDir;

        uniform mat4 model;
        uniform mat4 view;
        uniform mat4 projection;
        uniform vec3 viewPos;
        // Pre-computed normal matrix (transpose(inverse(model))) passed as uniform
        // to avoid per-vertex matrix inversion on GPU
        uniform mat3 normalMatrix;

        void main() {
            vec4 worldPos = model * vec4(aPos, 1.0);
            FragPos  = worldPos.xyz;
            Normal   = normalMatrix * aNormal;
            ViewDir  = normalize(viewPos - FragPos);
            gl_Position = projection * view * worldPos;
        }
    )";

    const char* accumFragSrc = R"(
        #version 330 core
        in vec3 FragPos;
        in vec3 Normal;
        in vec3 ViewDir;

        // MRT outputs:
        //   layout 0 -> accum texture  (RGBA16F)  — premultiplied-alpha weighted color
        //   layout 1 -> reveal texture (R16F)      — product of (1-alpha) weights
        layout (location = 0) out vec4 accumColor;
        layout (location = 1) out float revealage;

        uniform vec3  structureColor;
        uniform float opacity;
        uniform vec3  lightDir1;
        uniform vec3  lightDir2;

        void main() {
            vec3 norm = normalize(Normal);
            vec3 V    = normalize(ViewDir);

            // Two-sided: flip normal when back-facing
            if (dot(norm, V) < 0.0) norm = -norm;

            // Blinn-Phong lighting
            float ambientStrength = 0.7;
            vec3 ambient = ambientStrength * structureColor;

            vec3 L1    = normalize(lightDir1);
            float diff1 = max(dot(norm, L1), 0.0);
            vec3 H1    = normalize(L1 + V);
            float spec1 = pow(max(dot(norm, H1), 0.0), 32.0);

            vec3 L2    = normalize(lightDir2);
            float diff2 = max(dot(norm, L2), 0.0);
            vec3 H2    = normalize(L2 + V);
            float spec2 = pow(max(dot(norm, H2), 0.0), 16.0);

            vec3 diffuse  = (0.5 * diff1 + 0.3 * diff2) * structureColor;
            vec3 specular = (0.15 * spec1 + 0.08 * spec2) * vec3(1.0);

            float rim      = 1.0 - max(dot(V, norm), 0.0);
            rim            = smoothstep(0.3, 1.0, rim);
            vec3 rimColor  = rim * 0.2 * structureColor;

            vec3 color  = min(ambient + diffuse + specular + rimColor, vec3(1.0));
            float depthFade = clamp(1.0 - gl_FragCoord.z, 0.2, 1.0);
            float alpha = opacity * depthFade;

            // ── WBOIT weight function ──────────────────────────────────────
            // McGuire & Bavoil (2013) weight: balances near/far fragments.
            // Clamp depth to [0,1] from gl_FragCoord.z.
            float depth  = gl_FragCoord.z;
            // Weight emphasises closer surfaces and avoids near-zero for far ones.
            float weight = clamp(
                pow(alpha, 1.0) * 8.0 / (1e-3 + pow(depth / 200.0, 4.0)),
                1e-2, 3e3
            );

            // Output premultiplied-alpha color scaled by weight
            accumColor = vec4(color * alpha * weight, alpha * weight);
            revealage  = alpha; // compositor reads (1-alpha) product
        }
    )";

    // ── Composite pass shader ────────────────────────────────────────────────
    // Full-screen quad reads the two WBOIT textures and resolves the final color.
    const char* compositeVertSrc = R"(
        #version 330 core
        layout (location = 0) in vec2 aPos;
        out vec2 TexCoords;
        void main() {
            TexCoords   = aPos * 0.5 + 0.5;
            gl_Position = vec4(aPos, 0.0, 1.0);
        }
    )";

    const char* compositeFragSrc = R"(
        #version 330 core
        in vec2 TexCoords;
        out vec4 FragColor;

        uniform sampler2D accumTexture;
        uniform sampler2D revealTexture;

        void main() {
            vec4  accum    = texture(accumTexture,  TexCoords);
            float reveal   = texture(revealTexture, TexCoords).r;

            // Avoid divide-by-zero for fully transparent pixels
            if (accum.a < 1e-4) {
                discard;
            }

            // Reconstruct weighted-average color, then blend with background
            // using the revealage (accumulated transmittance).
            vec3 avgColor  = accum.rgb / accum.a;
            float transmit = clamp(1.0 - reveal, 0.0, 1.0);

            // Final over-compositing: color * (1 - transmittance) already baked in
            FragColor = vec4(avgColor, 1.0 - transmit);
        }
    )";

    // Helper lambda to compile a shader stage
    auto compileShader = [](GLenum type, const char* src) -> GLuint {
        GLuint shader = glCreateShader(type);
        glShaderSource(shader, 1, &src, nullptr);
        glCompileShader(shader);
        GLint ok = 0;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &ok);
        if (!ok) {
            char log[512];
            glGetShaderInfoLog(shader, 512, nullptr, log);
            Logger::error(std::string("Shader compile error: ") + log);
        }
        return shader;
    };

    // Build accumulation program
    {
        GLuint vs = compileShader(GL_VERTEX_SHADER,   accumVertSrc);
        GLuint fs = compileShader(GL_FRAGMENT_SHADER, accumFragSrc);
        m_accumProgram = glCreateProgram();
        glAttachShader(m_accumProgram, vs);
        glAttachShader(m_accumProgram, fs);
        glLinkProgram(m_accumProgram);
        glDeleteShader(vs);
        glDeleteShader(fs);

        // Cache uniform locations
        m_uModel         = glGetUniformLocation(m_accumProgram, "model");
        m_uView          = glGetUniformLocation(m_accumProgram, "view");
        m_uProjection    = glGetUniformLocation(m_accumProgram, "projection");
        m_uViewPos       = glGetUniformLocation(m_accumProgram, "viewPos");
        m_uNormalMatrix  = glGetUniformLocation(m_accumProgram, "normalMatrix");
        m_uLightDir1     = glGetUniformLocation(m_accumProgram, "lightDir1");
        m_uLightDir2     = glGetUniformLocation(m_accumProgram, "lightDir2");
        m_uColor         = glGetUniformLocation(m_accumProgram, "structureColor");
        m_uOpacity       = glGetUniformLocation(m_accumProgram, "opacity");
    }

    // Build composite program
    {
        GLuint vs = compileShader(GL_VERTEX_SHADER,   compositeVertSrc);
        GLuint fs = compileShader(GL_FRAGMENT_SHADER, compositeFragSrc);
        m_compositeProgram = glCreateProgram();
        glAttachShader(m_compositeProgram, vs);
        glAttachShader(m_compositeProgram, fs);
        glLinkProgram(m_compositeProgram);
        glDeleteShader(vs);
        glDeleteShader(fs);

        glUseProgram(m_compositeProgram);
        glUniform1i(glGetUniformLocation(m_compositeProgram, "accumTexture"),  0);
        glUniform1i(glGetUniformLocation(m_compositeProgram, "revealTexture"), 1);
        glUseProgram(0);
    }

    // Full-screen quad VBO (NDC, two triangles)
    static const float quadVerts[] = {
        -1.f, -1.f,   1.f, -1.f,   1.f,  1.f,
        -1.f, -1.f,   1.f,  1.f,  -1.f,  1.f
    };
    glGenVertexArrays(1, &m_quadVao);
    glGenBuffers(1, &m_quadVbo);
    glBindVertexArray(m_quadVao);
    glBindBuffer(GL_ARRAY_BUFFER, m_quadVbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVerts), quadVerts, GL_STATIC_DRAW);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    // WBOIT framebuffer will be created/resized lazily in render()
    m_wboitWidth  = 0;
    m_wboitHeight = 0;
}

// ──────────────────────────────────────────────────────────────────────────────
// WBOIT framebuffer helpers
// ──────────────────────────────────────────────────────────────────────────────

void StructureRenderer::ensureWboitFBO(int width, int height) {
    if (m_wboitFBO && m_wboitWidth == width && m_wboitHeight == height) return;

    // Destroy old resources
    destroyWboitFBO();

    m_wboitWidth  = width;
    m_wboitHeight = height;

    glGenFramebuffers(1, &m_wboitFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, m_wboitFBO);

    // Accumulation texture: RGBA16F (premultiplied weighted color + weight sum)
    glGenTextures(1, &m_accumTex);
    glBindTexture(GL_TEXTURE_2D, m_accumTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA16F, width, height, 0, GL_RGBA, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_accumTex, 0);

    // Revealage texture: R16F (transmittance accumulation)
    glGenTextures(1, &m_revealTex);
    glBindTexture(GL_TEXTURE_2D, m_revealTex);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R16F, width, height, 0, GL_RED, GL_FLOAT, nullptr);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, m_revealTex, 0);

    // Depth renderbuffer (shared with the main pass so structures occlude correctly)
    glGenRenderbuffers(1, &m_depthRBO);
    glBindRenderbuffer(GL_RENDERBUFFER, m_depthRBO);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, width, height);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_depthRBO);

    GLenum drawBuffers[2] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1 };
    glDrawBuffers(2, drawBuffers);

    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        Logger::error("WBOIT FBO incomplete!");

    glBindFramebuffer(GL_FRAMEBUFFER, 0);
}

void StructureRenderer::destroyWboitFBO() {
    if (m_wboitFBO)   { glDeleteFramebuffers(1,  &m_wboitFBO);   m_wboitFBO   = 0; }
    if (m_accumTex)   { glDeleteTextures(1,      &m_accumTex);   m_accumTex   = 0; }
    if (m_revealTex)  { glDeleteTextures(1,      &m_revealTex);  m_revealTex  = 0; }
    if (m_depthRBO)   { glDeleteRenderbuffers(1, &m_depthRBO);   m_depthRBO   = 0; }
}

// ──────────────────────────────────────────────────────────────────────────────

void StructureRenderer::setPatientData(PatientData* data) {
    if (m_patientData == data) return;
    m_patientData = data;
    m_needsRebuild = true;
}

// ──────────────────────────────────────────────────────────────────────────────
// Voxelization helpers
// ──────────────────────────────────────────────────────────────────────────────

// Point-in-polygon test (2D, ray casting)
static bool pointInPolygon2D(float px, float py,
                              const std::vector<std::array<double,3>>& poly,
                              int ax1, int ax2) {
    bool inside = false;
    size_t n = poly.size();
    for (size_t i = 0, j = n - 1; i < n; j = i++) {
        float yi = (float)poly[i][ax2], yj = (float)poly[j][ax2];
        float xi = (float)poly[i][ax1], xj = (float)poly[j][ax1];
        if (((yi > py) != (yj > py)) &&
            (px < (xj - xi) * (py - yi) / (yj - yi) + xi))
            inside = !inside;
    }
    return inside;
}

void StructureRenderer::tessellateStructure(size_t structureIndex) {
    auto* structures = m_patientData->getStructureSet();
    const auto* structure = structures->getStructure(structureIndex);
    if (!structure) {
        StructureMesh emptyMesh;
        emptyMesh.indexCount = 0;
        emptyMesh.center = glm::vec3(0.0f);
        emptyMesh.vao = emptyMesh.vbo = emptyMesh.ebo = 0;
        m_meshes.push_back(emptyMesh);
        return;
    }

    const auto& contours = structure->getContours();

    if (contours.size() < 2) {
        StructureMesh mesh;
        auto color = structure->getColor();
        mesh.color   = glm::vec3(color[0]/255.f, color[1]/255.f, color[2]/255.f);
        mesh.visible = structure->isVisible();
        mesh.indexCount = 0;
        mesh.center = glm::vec3(0.0f);
        mesh.vao = mesh.vbo = mesh.ebo = 0;
        m_meshes.push_back(mesh);
        Logger::info("    -> Skipped (need ≥2 contours for 3D mesh)");
        return;
    }

    // 1. Bounding box
    float minX = 1e30f, minY = 1e30f, minZ = 1e30f;
    float maxX = -1e30f, maxY = -1e30f, maxZ = -1e30f;
    for (const auto& c : contours) {
        for (const auto& pt : c.points) {
            minX = std::min(minX,(float)pt[0]); maxX = std::max(maxX,(float)pt[0]);
            minY = std::min(minY,(float)pt[1]); maxY = std::max(maxY,(float)pt[1]);
            minZ = std::min(minZ,(float)pt[2]); maxZ = std::max(maxZ,(float)pt[2]);
        }
    }
    float margin = 2.0f;
    minX -= margin; minY -= margin; minZ -= margin;
    maxX += margin; maxY += margin; maxZ += margin;

    // 2. Voxel size — adaptive
    float rangeX = maxX-minX, rangeY = maxY-minY, rangeZ = maxZ-minZ;
    float maxRange = std::max({rangeX, rangeY, rangeZ});
    float voxelSize = maxRange / 100.0f;
    if (voxelSize < 0.5f) voxelSize = 0.5f;

    int nx = std::max(2,(int)std::ceil(rangeX/voxelSize)+1);
    int ny = std::max(2,(int)std::ceil(rangeY/voxelSize)+1);
    int nz = std::max(2,(int)std::ceil(rangeZ/voxelSize)+1);

    if ((long long)nx*ny*nz > 8000000LL) {
        voxelSize *= std::cbrt((double)nx*ny*nz/8000000.0);
        nx = std::max(2,(int)std::ceil(rangeX/voxelSize)+1);
        ny = std::max(2,(int)std::ceil(rangeY/voxelSize)+1);
        nz = std::max(2,(int)std::ceil(rangeZ/voxelSize)+1);
    }

    // 3. Group contours by Z slice
    std::unordered_map<int, std::vector<size_t>> sliceContours;
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        if (contours[ci].points.empty()) continue;
        float z  = contours[ci].points[0][2];
        int  iz  = (int)std::round((z - minZ) / voxelSize);
        sliceContours[iz].push_back(ci);
    }

    // 4. Voxelize
    std::vector<float> volume(nx*ny*nz, 0.0f);
    auto idx = [&](int x,int y,int z) -> size_t {
        return (size_t)z*ny*nx + (size_t)y*nx + x;
    };

    for (auto& [iz, cIndices] : sliceContours) {
        if (iz < 0 || iz >= nz) continue;
        for (int iy = 0; iy < ny; ++iy) {
            float py = minY + iy*voxelSize;
            for (int ix = 0; ix < nx; ++ix) {
                float px = minX + ix*voxelSize;
                for (size_t ci : cIndices) {
                    if (contours[ci].points.size() < 3) continue;
                    if (pointInPolygon2D(px, py, contours[ci].points, 0, 1)) {
                        volume[idx(ix,iy,iz)] = 1.0f;
                        break;
                    }
                }
            }
        }
    }

    // 4b. Interpolate between slices
    std::vector<int> sortedZ;
    sortedZ.reserve(sliceContours.size());
    for (auto& [iz, _] : sliceContours) sortedZ.push_back(iz);
    std::sort(sortedZ.begin(), sortedZ.end());

    for (size_t si = 0; si+1 < sortedZ.size(); ++si) {
        int z0 = sortedZ[si], z1 = sortedZ[si+1];
        if (z1-z0 <= 1) continue;
        for (int iz = z0+1; iz < z1; ++iz)
            for (int iy = 0; iy < ny; ++iy)
                for (int ix = 0; ix < nx; ++ix)
                    if (volume[idx(ix,iy,z0)] > 0.5f && volume[idx(ix,iy,z1)] > 0.5f)
                        volume[idx(ix,iy,iz)] = 1.0f;
    }

    // 5. Gaussian blur
    {
        std::vector<float> blurred(nx*ny*nz, 0.0f);
        for (int iz = 1; iz < nz-1; ++iz)
            for (int iy = 1; iy < ny-1; ++iy)
                for (int ix = 1; ix < nx-1; ++ix) {
                    float sum=0.f, wTotal=0.f;
                    for (int dz=-1;dz<=1;++dz)
                        for (int dy=-1;dy<=1;++dy)
                            for (int dx=-1;dx<=1;++dx) {
                                float w = 1.f/(1.f+std::abs(dx)+std::abs(dy)+std::abs(dz));
                                sum    += w*volume[idx(ix+dx,iy+dy,iz+dz)];
                                wTotal += w;
                            }
                    blurred[idx(ix,iy,iz)] = sum/wTotal;
                }
        volume = std::move(blurred);
    }

    // 6. Marching Cubes at isovalue 0.5
    float isovalue = 0.5f;
    float scale    = 0.001f; // mm -> world units

    std::vector<float>        vertices;
    std::vector<unsigned int> indices;

    auto getVal = [&](int x,int y,int z) -> float {
        if (x<0||x>=nx||y<0||y>=ny||z<0||z>=nz) return 0.0f;
        return volume[idx(x,y,z)];
    };

    auto computeGradient = [&](float fx,float fy,float fz) -> glm::vec3 {
        float eps = 0.5f;
        auto sample = [&](float sx,float sy,float sz) -> float {
            int ix=(int)sx, iy=(int)sy, iz=(int)sz;
            float lx=sx-ix, ly=sy-iy, lz=sz-iz;
            auto s=[&](int x,int y,int z){ return getVal(x,y,z); };
            float c000=s(ix,iy,iz),   c100=s(ix+1,iy,iz);
            float c010=s(ix,iy+1,iz), c110=s(ix+1,iy+1,iz);
            float c001=s(ix,iy,iz+1), c101=s(ix+1,iy,iz+1);
            float c011=s(ix,iy+1,iz+1),c111=s(ix+1,iy+1,iz+1);
            float c00=c000*(1-lx)+c100*lx, c10=c010*(1-lx)+c110*lx;
            float c01=c001*(1-lx)+c101*lx, c11=c011*(1-lx)+c111*lx;
            float c0=c00*(1-ly)+c10*ly,    c1=c01*(1-ly)+c11*ly;
            return c0*(1-lz)+c1*lz;
        };
        float gx=sample(fx+eps,fy,fz)-sample(fx-eps,fy,fz);
        float gy=sample(fx,fy+eps,fz)-sample(fx,fy-eps,fz);
        float gz=sample(fx,fy,fz+eps)-sample(fx,fy,fz-eps);
        glm::vec3 g(-gx,-gy,-gz);
        float len=glm::length(g);
        return len>1e-8f ? g/len : glm::vec3(0,0,1);
    };

    // Edge vertex deduplication.
    // FIX: The original code used (ix, iy, iz, e) as key, which is NOT unique
    // across cells sharing an edge — adjacent cells would create duplicate vertices
    // at the same position with different normals, causing lighting flicker.
    //
    // Correct approach: canonicalize each edge by the absolute grid positions of
    // its two endpoints (always ordered so lower coord comes first). This guarantees
    // exactly one vertex per edge regardless of which cell processes it first.
    struct EdgeKey {
        int ax0,ay0,az0, ax1,ay1,az1;
        bool operator==(const EdgeKey& o) const {
            return ax0==o.ax0&&ay0==o.ay0&&az0==o.az0&&
                   ax1==o.ax1&&ay1==o.ay1&&az1==o.az1;
        }
    };
    struct EdgeKeyHash {
        size_t operator()(const EdgeKey& k) const {
            // FNV-style mix of all 6 ints
            size_t h = 2166136261u;
            auto mix = [&](int v){ h = (h^(size_t)v)*16777619u; };
            mix(k.ax0); mix(k.ay0); mix(k.az0);
            mix(k.ax1); mix(k.ay1); mix(k.az1);
            return h;
        }
    };
    std::unordered_map<EdgeKey, unsigned int, EdgeKeyHash> edgeVertexMap;

    for (int iz = 0; iz < nz-1; ++iz) {
        for (int iy = 0; iy < ny-1; ++iy) {
            for (int ix = 0; ix < nx-1; ++ix) {
                float cornerVals[8];
                cornerVals[0]=getVal(ix,  iy,  iz);
                cornerVals[1]=getVal(ix+1,iy,  iz);
                cornerVals[2]=getVal(ix+1,iy+1,iz);
                cornerVals[3]=getVal(ix,  iy+1,iz);
                cornerVals[4]=getVal(ix,  iy,  iz+1);
                cornerVals[5]=getVal(ix+1,iy,  iz+1);
                cornerVals[6]=getVal(ix+1,iy+1,iz+1);
                cornerVals[7]=getVal(ix,  iy+1,iz+1);

                int cubeIndex = 0;
                for (int c=0;c<8;++c)
                    if (cornerVals[c]>=isovalue) cubeIndex|=(1<<c);

                if (mc::edgeTable[cubeIndex]==0) continue;

                unsigned int edgeVerts[12]={};
                for (int e=0;e<12;++e) {
                    if (!(mc::edgeTable[cubeIndex]&(1<<e))) continue;

                    int c0=mc::edgeConnection[e][0];
                    int c1=mc::edgeConnection[e][1];

                    // Absolute grid positions of the two edge endpoints
                    int ax0=ix+mc::cornerOffsets[c0][0];
                    int ay0=iy+mc::cornerOffsets[c0][1];
                    int az0=iz+mc::cornerOffsets[c0][2];
                    int ax1=ix+mc::cornerOffsets[c1][0];
                    int ay1=iy+mc::cornerOffsets[c1][1];
                    int az1=iz+mc::cornerOffsets[c1][2];

                    // Canonical ordering: lower endpoint first (lexicographic)
                    if (std::tie(ax0,ay0,az0)>std::tie(ax1,ay1,az1)) {
                        std::swap(ax0,ax1);
                        std::swap(ay0,ay1);
                        std::swap(az0,az1);
                        std::swap(c0,c1); // keep c0/c1 consistent for interpolation dir
                    }

                    EdgeKey key{ax0,ay0,az0,ax1,ay1,az1};
                    auto it=edgeVertexMap.find(key);
                    if (it!=edgeVertexMap.end()) {
                        edgeVerts[e]=it->second;
                    } else {
                        float v0=cornerVals[c0], v1=cornerVals[c1];
                        float t=(std::abs(v1-v0)>1e-10f)?(isovalue-v0)/(v1-v0):0.5f;
                        t=std::clamp(t,0.0f,1.0f);

                        float gx_=ix+mc::cornerOffsets[c0][0]+t*(mc::cornerOffsets[c1][0]-mc::cornerOffsets[c0][0]);
                        float gy_=iy+mc::cornerOffsets[c0][1]+t*(mc::cornerOffsets[c1][1]-mc::cornerOffsets[c0][1]);
                        float gz_=iz+mc::cornerOffsets[c0][2]+t*(mc::cornerOffsets[c1][2]-mc::cornerOffsets[c0][2]);

                        float wx=(minX+gx_*voxelSize)*scale;
                        float wy=(minY+gy_*voxelSize)*scale;
                        float wz=(minZ+gz_*voxelSize)*scale;

                        glm::vec3 n=computeGradient(gx_,gy_,gz_);

                        unsigned int vi=(unsigned int)(vertices.size()/6);
                        vertices.push_back(wx); vertices.push_back(wy); vertices.push_back(wz);
                        vertices.push_back(n.x); vertices.push_back(n.y); vertices.push_back(n.z);

                        edgeVertexMap[key]=vi;
                        edgeVerts[e]=vi;
                    }
                }

                for (int t=0;mc::triTable[cubeIndex][t]!=-1;t+=3) {
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t]]);
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t+1]]);
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t+2]]);
                }
            }
        }
    }

    if (vertices.empty()||indices.empty()) {
        StructureMesh mesh;
        auto color=structure->getColor();
        mesh.color=glm::vec3(color[0]/255.f,color[1]/255.f,color[2]/255.f);
        mesh.visible=structure->isVisible();
        mesh.indexCount=0; mesh.center=glm::vec3(0.f);
        mesh.vao=mesh.vbo=mesh.ebo=0;
        m_meshes.push_back(mesh);
        Logger::info("    -> Skipped (marching cubes produced no geometry)");
        return;
    }

    // 7. Smooth normals by averaging face normals per vertex
    {
        size_t vertCount=vertices.size()/6;
        std::vector<glm::vec3> smoothNormals(vertCount,glm::vec3(0.f));

        for (size_t i=0;i+2<indices.size();i+=3) {
            unsigned int i0=indices[i],i1=indices[i+1],i2=indices[i+2];
            glm::vec3 v0(vertices[i0*6],vertices[i0*6+1],vertices[i0*6+2]);
            glm::vec3 v1(vertices[i1*6],vertices[i1*6+1],vertices[i1*6+2]);
            glm::vec3 v2(vertices[i2*6],vertices[i2*6+1],vertices[i2*6+2]);
            glm::vec3 faceN=glm::cross(v1-v0,v2-v0);
            smoothNormals[i0]+=faceN;
            smoothNormals[i1]+=faceN;
            smoothNormals[i2]+=faceN;
        }
        for (size_t i=0;i<vertCount;++i) {
            glm::vec3 n=glm::normalize(smoothNormals[i]);
            if (glm::length(smoothNormals[i])<1e-8f) n=glm::vec3(0,0,1);
            vertices[i*6+3]=n.x; vertices[i*6+4]=n.y; vertices[i*6+5]=n.z;
        }
    }

    // Bounding box center
    glm::vec3 minBound(1e30f),maxBound(-1e30f);
    size_t vertCount=vertices.size()/6;
    for (size_t i=0;i<vertCount;++i) {
        glm::vec3 pos(vertices[i*6],vertices[i*6+1],vertices[i*6+2]);
        minBound=glm::min(minBound,pos);
        maxBound=glm::max(maxBound,pos);
    }

    // 8. Upload to GPU
    StructureMesh mesh;
    auto color=structure->getColor();
    mesh.color=glm::vec3(color[0]/255.f,color[1]/255.f,color[2]/255.f);
    mesh.visible=structure->isVisible();
    mesh.indexCount=indices.size();
    mesh.center=(minBound+maxBound)*0.5f;

    glGenVertexArrays(1,&mesh.vao);
    glGenBuffers(1,&mesh.vbo);
    glGenBuffers(1,&mesh.ebo);
    glBindVertexArray(mesh.vao);
    glBindBuffer(GL_ARRAY_BUFFER,mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER,vertices.size()*sizeof(float),vertices.data(),GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER,indices.size()*sizeof(unsigned int),indices.data(),GL_STATIC_DRAW);
    glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,6*sizeof(float),(void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1,3,GL_FLOAT,GL_FALSE,6*sizeof(float),(void*)(3*sizeof(float)));
    glEnableVertexAttribArray(1);
    glBindVertexArray(0);

    m_meshes.push_back(mesh);
}

// ──────────────────────────────────────────────────────────────────────────────

void StructureRenderer::buildMeshes() {
    std::lock_guard<std::mutex> lock(m_meshMutex);
    buildMeshes_unlocked();
}

void StructureRenderer::buildMeshes_unlocked() {
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1,&mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1,&mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1,&mesh.ebo);
    }
    m_meshes.clear();

    if (!m_patientData||!m_patientData->getStructureSet()) return;

    auto* structures=m_patientData->getStructureSet();
    size_t structureCount=structures->getCount();
    Logger::info("Building meshes for "+std::to_string(structureCount)+" structures");

    for (size_t i=0;i<structureCount;++i) {
        const auto* structure=structures->getStructure(i);
        if (structure)
            Logger::info("  Structure "+std::to_string(i)+": "+structure->getName()+
                         " ("+std::to_string(structure->getContourCount())+" contours)");
        tessellateStructure(i);
    }
    Logger::info("Built "+std::to_string(m_meshes.size())+" meshes");

    if (m_meshes.size()!=structureCount)
        Logger::error("MISMATCH: "+std::to_string(m_meshes.size())+" meshes vs "+
                      std::to_string(structureCount)+" structures!");

    m_needsRebuild=false;
}

// ──────────────────────────────────────────────────────────────────────────────
// render() — Weighted Blended Order-Independent Transparency (WBOIT)
//
// Algorithm (McGuire & Bavoil, 2013):
//   Pass 1 — Accumulation:
//       Render all transparent meshes into two MRT textures:
//         accumTex  (RGBA16F): sum of (color*alpha*w, alpha*w)  [additive blend]
//         revealTex (R16F):    product of (1-alpha) per pixel   [multiplicative blend]
//       No depth writes, no sorting required.
//   Pass 2 — Composite:
//       Full-screen quad reads both textures and resolves final color onto the
//       default framebuffer using normal over-blending.
// ──────────────────────────────────────────────────────────────────────────────

void StructureRenderer::render(const glm::mat4& view, const glm::mat4& projection) {
    std::lock_guard<std::mutex> lock(m_meshMutex);

    if (!m_patientData||!m_patientData->getStructureSet()) return;
    if (m_needsRebuild) buildMeshes_unlocked();
    if (m_meshes.empty()) return;

    auto* structures=m_patientData->getStructureSet();
    size_t structureCount=structures->getCount();

    if (m_meshes.size()!=structureCount) {
        Logger::warn("Mesh/structure count mismatch during render, rebuilding...");
        buildMeshes_unlocked();
        if (m_meshes.size()!=structureCount) {
            Logger::error("Failed to fix mesh/structure alignment!");
            return;
        }
    }

    // Query current viewport to size the WBOIT FBO
    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT, viewport);
    int vpW=viewport[2], vpH=viewport[3];
    ensureWboitFBO(vpW, vpH);

    // -90° rotation to align with anatomical axes
    glm::mat4 model = glm::rotate(glm::mat4(1.f), -glm::half_pi<float>(), glm::vec3(1,0,0));

    // Normal matrix (pre-computed once per frame, uploaded as uniform — avoids
    // per-vertex matrix inversion in the vertex shader)
    glm::mat3 normalMatrix = glm::mat3(glm::transpose(glm::inverse(model)));

    // Camera position
    glm::mat4 invView = glm::inverse(view);
    glm::vec3 viewPos(invView[3][0], invView[3][1], invView[3][2]);

    // Lights (world-space, fixed directions)
    glm::vec3 lightDir1 = glm::normalize(glm::vec3( 0.5f,  0.8f, 0.6f));
    glm::vec3 lightDir2 = glm::normalize(glm::vec3(-0.4f, -0.3f, 0.8f));

    // ── Copy opaque depth into the WBOIT FBO ──────────────────────────────
    // Blit the current depth buffer so that structures are correctly occluded
    // by opaque geometry rendered before this call.
    GLint readFBO=0;
    glGetIntegerv(GL_READ_FRAMEBUFFER_BINDING, &readFBO);
    glBindFramebuffer(GL_READ_FRAMEBUFFER, readFBO);
    glBindFramebuffer(GL_DRAW_FRAMEBUFFER, m_wboitFBO);
    glBlitFramebuffer(0,0,vpW,vpH, 0,0,vpW,vpH, GL_DEPTH_BUFFER_BIT, GL_NEAREST);

    // ── Pass 1: Accumulation ───────────────────────────────────────────────
    glBindFramebuffer(GL_FRAMEBUFFER, m_wboitFBO);

    // Clear accum to (0,0,0,0) and revealage to (1) — 1 = no occlusion yet
    static const float zeroClear[4] = {0,0,0,0};
    static const float oneClear[4]  = {1,1,1,1};
    glClearBufferfv(GL_COLOR, 0, zeroClear); // accumTex
    glClearBufferfv(GL_COLOR, 1, oneClear);  // revealTex

    // Depth test ON (respect opaque geometry), depth write OFF (all layers contribute)
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    glDisable(GL_CULL_FACE); // render both faces without two-pass trick;
                              // two-sided lighting is handled in the shader

    glEnable(GL_BLEND);
    // accumTex  attachment: GL_ONE / GL_ONE  (additive accumulation)
    glBlendFunci(0, GL_ONE, GL_ONE);
    // revealTex attachment: GL_ZERO / GL_ONE_MINUS_SRC_COLOR  (multiplicative)
    glBlendFunci(1, GL_ZERO, GL_ONE_MINUS_SRC_COLOR);

    glUseProgram(m_accumProgram);

    // Upload per-frame uniforms (cached locations — no glGetUniformLocation each frame)
    glUniformMatrix4fv(m_uModel,        1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(m_uView,         1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(m_uProjection,   1, GL_FALSE, glm::value_ptr(projection));
    glUniform3fv      (m_uViewPos,      1, glm::value_ptr(viewPos));
    glUniformMatrix3fv(m_uNormalMatrix, 1, GL_FALSE, glm::value_ptr(normalMatrix));
    glUniform3fv      (m_uLightDir1,    1, glm::value_ptr(lightDir1));
    glUniform3fv      (m_uLightDir2,    1, glm::value_ptr(lightDir2));

    // Draw all visible meshes — no sorting needed with WBOIT
    for (size_t i = 0; i < m_meshes.size(); ++i) {
        const auto& mesh = m_meshes[i];
        if (mesh.indexCount == 0) continue;
        const auto* structure = structures->getStructure(i);
        float opacity = 0.03f;

        if (structure->getName().find("GTV") != std::string::npos)
            opacity = 0.30f; // target
        else if (structure->getName().find("CTV") != std::string::npos)
            opacity = 0.20f; // clinical target
        else if (structure->getName().find("PTV") != std::string::npos)
            opacity = 0.15; // high-risk
        else if (structure->getName().find("BODY") != std::string::npos)
            opacity = 0.01f; // outer shell
        else
                opacity = 0.03f; // default for other structures

        glUniform1f(m_uOpacity, opacity);
        if (!structure || !structure->isVisible()) continue;

        glUniform3fv(m_uColor, 1, glm::value_ptr(mesh.color));
        glBindVertexArray(mesh.vao);
        glDrawElements(GL_TRIANGLES, (GLsizei)mesh.indexCount, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
    glUseProgram(0);

    // ── Pass 2: Composite onto default framebuffer ─────────────────────────
    glBindFramebuffer(GL_FRAMEBUFFER, (GLuint)readFBO);

    // Standard over-blending for the resolved transparent layer
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendEquation(GL_FUNC_ADD);

    // Depth test OFF — the composite quad covers the full screen
    glDisable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);

    glUseProgram(m_compositeProgram);
    glActiveTexture(GL_TEXTURE0); glBindTexture(GL_TEXTURE_2D, m_accumTex);
    glActiveTexture(GL_TEXTURE1); glBindTexture(GL_TEXTURE_2D, m_revealTex);

    glBindVertexArray(m_quadVao);
    glDrawArrays(GL_TRIANGLES, 0, 6);
    glBindVertexArray(0);

    glUseProgram(0);

    // ── Restore GL state ───────────────────────────────────────────────────
    glDepthMask(GL_TRUE);
    glEnable(GL_DEPTH_TEST);
    glDisable(GL_BLEND);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);
    glActiveTexture(GL_TEXTURE0);
}

// ──────────────────────────────────────────────────────────────────────────────

void StructureRenderer::cleanup() {
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1,&mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1,&mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1,&mesh.ebo);
    }
    m_meshes.clear();

    destroyWboitFBO();

    if (m_quadVao) { glDeleteVertexArrays(1,&m_quadVao); m_quadVao=0; }
    if (m_quadVbo) { glDeleteBuffers(1,&m_quadVbo);      m_quadVbo=0; }

    if (m_accumProgram)     { glDeleteProgram(m_accumProgram);     m_accumProgram=0;     }
    if (m_compositeProgram) { glDeleteProgram(m_compositeProgram); m_compositeProgram=0; }
}

} // namespace optirad