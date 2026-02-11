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
    // Blinn-Phong shader with very bright, saturated colors
    const char* vertexShaderSource = R"(
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

        void main() {
            vec4 worldPos = model * vec4(aPos, 1.0);
            FragPos = worldPos.xyz;
            Normal = mat3(transpose(inverse(model))) * aNormal;
            ViewDir = normalize(viewPos - FragPos);
            gl_Position = projection * view * worldPos;
        }
    )";

    const char* fragmentShaderSource = R"(
        #version 330 core
        in vec3 FragPos;
        in vec3 Normal;
        in vec3 ViewDir;
        out vec4 FragColor;

        uniform vec3 structureColor;
        uniform float opacity;
        uniform vec3 lightDir1;
        uniform vec3 lightDir2;

        void main() {
            vec3 norm = normalize(Normal);
            // Two-sided lighting: flip normal if facing away from camera
            if (dot(norm, ViewDir) < 0.0)
                norm = -norm;

            vec3 V = normalize(ViewDir);

            // Very high ambient for bright base color
            float ambientStrength = 0.7;
            vec3 ambient = ambientStrength * structureColor;

            // Bright diffuse lighting
            vec3 L1 = normalize(lightDir1);
            float diff1 = max(dot(norm, L1), 0.0);
            vec3 H1 = normalize(L1 + V);
            float spec1 = pow(max(dot(norm, H1), 0.0), 32.0);

            vec3 L2 = normalize(lightDir2);
            float diff2 = max(dot(norm, L2), 0.0);
            vec3 H2 = normalize(L2 + V);
            float spec2 = pow(max(dot(norm, H2), 0.0), 16.0);

            // Strong diffuse contribution
            vec3 diffuse = (0.5 * diff1 + 0.3 * diff2) * structureColor;
            
            // Subtle white specular highlights
            vec3 specular = (0.15 * spec1 + 0.08 * spec2) * vec3(1.0);

            // Bright rim light for silhouettes
            float rim = 1.0 - max(dot(V, norm), 0.0);
            rim = smoothstep(0.3, 1.0, rim);
            vec3 rimColor = rim * 0.2 * structureColor;

            // Combine everything - NO tone mapping for maximum brightness
            vec3 result = ambient + diffuse + specular + rimColor;
            
            // Clamp to prevent over-saturation but keep it bright
            result = min(result, vec3(1.0));

            FragColor = vec4(result, opacity);
        }
    )";

    GLuint vs = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vs, 1, &vertexShaderSource, nullptr);
    glCompileShader(vs);

    GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fs, 1, &fragmentShaderSource, nullptr);
    glCompileShader(fs);

    m_shaderProgram = glCreateProgram();
    glAttachShader(m_shaderProgram, vs);
    glAttachShader(m_shaderProgram, fs);
    glLinkProgram(m_shaderProgram);
    glDeleteShader(vs);
    glDeleteShader(fs);
}

void StructureRenderer::setPatientData(PatientData* data) {
    if (m_patientData == data) return;
    m_patientData = data;
    m_needsRebuild = true;
}

// ──────────────────────── Voxelization helpers ────────────────────────

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

// Voxelize a single structure into a 3D binary grid, then run Marching Cubes
void StructureRenderer::tessellateStructure(size_t structureIndex) {
    auto* structures = m_patientData->getStructureSet();
    const auto* structure = structures->getStructure(structureIndex);
    if (!structure) {
        // Even if structure is null, add empty mesh to maintain index alignment (FIX ?)
        StructureMesh emptyMesh;
        emptyMesh.indexCount = 0;
        emptyMesh.center = glm::vec3(0.0f);
        emptyMesh.vao = 0;
        emptyMesh.vbo = 0;
        emptyMesh.ebo = 0;
        m_meshes.push_back(emptyMesh);
        return;
    }

    const auto& contours = structure->getContours();
    
    // Create empty mesh placeholder for structures with insufficient contours
    // This keeps mesh indices aligned with structure indices
    // Need at least 2 contours to interpolate between slices for a 3D volume
    if (contours.size() < 2) {
        StructureMesh mesh;
        auto color = structure->getColor();
        mesh.color = glm::vec3(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
        mesh.visible = structure->isVisible();
        mesh.indexCount = 0;
        mesh.center = glm::vec3(0.0f);
        mesh.vao = 0;
        mesh.vbo = 0;
        mesh.ebo = 0;
        m_meshes.push_back(mesh);
        Logger::info("    -> Skipped (need ≥2 contours for 3D mesh)");
        return;
    }

    // 1. Compute bounding box of all contour points (in mm)
    float minX = 1e30f, minY = 1e30f, minZ = 1e30f;
    float maxX = -1e30f, maxY = -1e30f, maxZ = -1e30f;
    for (const auto& c : contours) {
        for (const auto& pt : c.points) {
            minX = std::min(minX, (float)pt[0]); maxX = std::max(maxX, (float)pt[0]);
            minY = std::min(minY, (float)pt[1]); maxY = std::max(maxY, (float)pt[1]);
            minZ = std::min(minZ, (float)pt[2]); maxZ = std::max(maxZ, (float)pt[2]);
        }
    }

    // Add small margin
    float margin = 2.0f; // mm
    minX -= margin; minY -= margin; minZ -= margin;
    maxX += margin; maxY += margin; maxZ += margin;

    // 2. Choose voxel size — adaptive, aim for ~80-120 voxels along longest axis
    float rangeX = maxX - minX, rangeY = maxY - minY, rangeZ = maxZ - minZ;
    float maxRange = std::max({rangeX, rangeY, rangeZ});
    float voxelSize = maxRange / 100.0f;
    if (voxelSize < 0.5f) voxelSize = 0.5f;

    int nx = std::max(2, (int)std::ceil(rangeX / voxelSize) + 1);
    int ny = std::max(2, (int)std::ceil(rangeY / voxelSize) + 1);
    int nz = std::max(2, (int)std::ceil(rangeZ / voxelSize) + 1);

    // Cap grid size to prevent huge memory usage
    if ((long long)nx * ny * nz > 8000000LL) {
        voxelSize *= std::cbrt((double)nx * ny * nz / 8000000.0);
        nx = std::max(2, (int)std::ceil(rangeX / voxelSize) + 1);
        ny = std::max(2, (int)std::ceil(rangeY / voxelSize) + 1);
        nz = std::max(2, (int)std::ceil(rangeZ / voxelSize) + 1);
    }

    // 3. Group contours by Z slice
    std::unordered_map<int, std::vector<size_t>> sliceContours;
    for (size_t ci = 0; ci < contours.size(); ++ci) {
        if (contours[ci].points.empty()) continue;
        float z = contours[ci].points[0][2];
        int iz = (int)std::round((z - minZ) / voxelSize);
        sliceContours[iz].push_back(ci);
    }

    // 4. Voxelize: for each Z-slice that has contours, fill interior via point-in-polygon
    std::vector<float> volume(nx * ny * nz, 0.0f);

    auto idx = [&](int x, int y, int z) -> size_t {
        return (size_t)z * ny * nx + (size_t)y * nx + x;
    };

    for (auto& [iz, cIndices] : sliceContours) {
        if (iz < 0 || iz >= nz) continue;
        for (int iy = 0; iy < ny; ++iy) {
            float py = minY + iy * voxelSize;
            for (int ix = 0; ix < nx; ++ix) {
                float px = minX + ix * voxelSize;
                for (size_t ci : cIndices) {
                    if (contours[ci].points.size() < 3) continue;
                    if (pointInPolygon2D(px, py, contours[ci].points, 0, 1)) {
                        volume[idx(ix, iy, iz)] = 1.0f;
                        break; // inside at least one contour on this slice
                    }
                }
            }
        }
    }

    // 4b. Interpolate between slices: fill intermediate Z slices
    // Sort unique Z indices
    std::vector<int> sortedZ;
    sortedZ.reserve(sliceContours.size());
    for (auto& [iz, _] : sliceContours) sortedZ.push_back(iz);
    std::sort(sortedZ.begin(), sortedZ.end());

    for (size_t si = 0; si + 1 < sortedZ.size(); ++si) {
        int z0 = sortedZ[si], z1 = sortedZ[si + 1];
        if (z1 - z0 <= 1) continue;
        // Fill intermediate slices by interpolating (if both neighbors have a filled voxel, fill it)
        for (int iz = z0 + 1; iz < z1; ++iz) {
            for (int iy = 0; iy < ny; ++iy) {
                for (int ix = 0; ix < nx; ++ix) {
                    float v0 = volume[idx(ix, iy, z0)];
                    float v1 = volume[idx(ix, iy, z1)];
                    if (v0 > 0.5f && v1 > 0.5f) {
                        volume[idx(ix, iy, iz)] = 1.0f;
                    }
                }
            }
        }
    }

    // 5. Gaussian blur for smooth isosurface (3x3x3 kernel, single pass)
    {
        std::vector<float> blurred(nx * ny * nz, 0.0f);
        for (int iz = 1; iz < nz - 1; ++iz) {
            for (int iy = 1; iy < ny - 1; ++iy) {
                for (int ix = 1; ix < nx - 1; ++ix) {
                    float sum = 0.0f;
                    float wTotal = 0.0f;
                    for (int dz = -1; dz <= 1; ++dz) {
                        for (int dy = -1; dy <= 1; ++dy) {
                            for (int dx = -1; dx <= 1; ++dx) {
                                float w = 1.0f / (1.0f + std::abs(dx) + std::abs(dy) + std::abs(dz));
                                sum += w * volume[idx(ix+dx, iy+dy, iz+dz)];
                                wTotal += w;
                            }
                        }
                    }
                    blurred[idx(ix, iy, iz)] = sum / wTotal;
                }
            }
        }
        volume = std::move(blurred);
    }

    // 6. Marching Cubes at isovalue 0.5
    float isovalue = 0.5f;
    float scale = 0.001f; // mm -> world units

    std::vector<float> vertices;   // pos(3) + normal(3) per vertex
    std::vector<unsigned int> indices;

    // Hash for vertex deduplication
    struct Vec3Hash {
        size_t operator()(const std::tuple<int,int,int,int>& k) const {
            auto h = std::hash<int>();
            size_t seed = h(std::get<0>(k));
            seed ^= h(std::get<1>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            seed ^= h(std::get<2>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            seed ^= h(std::get<3>(k)) + 0x9e3779b9 + (seed<<6) + (seed>>2);
            return seed;
        }
    };
    // Map from (grid cell + edge index) -> vertex index
    std::unordered_map<std::tuple<int,int,int,int>, unsigned int, Vec3Hash> edgeVertexMap;

    auto getVal = [&](int x, int y, int z) -> float {
        if (x < 0 || x >= nx || y < 0 || y >= ny || z < 0 || z >= nz) return 0.0f;
        return volume[idx(x, y, z)];
    };

    // Compute gradient (central differences) for normals
    auto computeGradient = [&](float fx, float fy, float fz) -> glm::vec3 {
        // fx, fy, fz in grid-space (can be fractional)
        float eps = 0.5f;
        auto sample = [&](float sx, float sy, float sz) -> float {
            int ix = (int)sx, iy = (int)sy, iz = (int)sz;
            float lx = sx - ix, ly = sy - iy, lz = sz - iz;
            // Trilinear interpolation
            auto s = [&](int x, int y, int z) { return getVal(x,y,z); };
            float c000 = s(ix,iy,iz), c100 = s(ix+1,iy,iz);
            float c010 = s(ix,iy+1,iz), c110 = s(ix+1,iy+1,iz);
            float c001 = s(ix,iy,iz+1), c101 = s(ix+1,iy,iz+1);
            float c011 = s(ix,iy+1,iz+1), c111 = s(ix+1,iy+1,iz+1);
            float c00 = c000*(1-lx)+c100*lx;
            float c10 = c010*(1-lx)+c110*lx;
            float c01 = c001*(1-lx)+c101*lx;
            float c11 = c011*(1-lx)+c111*lx;
            float c0 = c00*(1-ly)+c10*ly;
            float c1 = c01*(1-ly)+c11*ly;
            return c0*(1-lz)+c1*lz;
        };
        float gx = sample(fx+eps,fy,fz) - sample(fx-eps,fy,fz);
        float gy = sample(fx,fy+eps,fz) - sample(fx,fy-eps,fz);
        float gz = sample(fx,fy,fz+eps) - sample(fx,fy,fz-eps);
        glm::vec3 g(-gx, -gy, -gz); // Negative gradient points outward
        float len = glm::length(g);
        return len > 1e-8f ? g / len : glm::vec3(0,0,1);
    };

    for (int iz = 0; iz < nz - 1; ++iz) {
        for (int iy = 0; iy < ny - 1; ++iy) {
            for (int ix = 0; ix < nx - 1; ++ix) {
                // Get values at 8 corners
                float cornerVals[8];
                cornerVals[0] = getVal(ix,   iy,   iz);
                cornerVals[1] = getVal(ix+1, iy,   iz);
                cornerVals[2] = getVal(ix+1, iy+1, iz);
                cornerVals[3] = getVal(ix,   iy+1, iz);
                cornerVals[4] = getVal(ix,   iy,   iz+1);
                cornerVals[5] = getVal(ix+1, iy,   iz+1);
                cornerVals[6] = getVal(ix+1, iy+1, iz+1);
                cornerVals[7] = getVal(ix,   iy+1, iz+1);

                int cubeIndex = 0;
                for (int c = 0; c < 8; ++c)
                    if (cornerVals[c] >= isovalue) cubeIndex |= (1 << c);

                // Validate marching cubes index
                if (cubeIndex < 0 || cubeIndex >= 256) {
                    Logger::error("Invalid marching cubes index: " + std::to_string(cubeIndex));
                    continue;
                }

                if (mc::edgeTable[cubeIndex] == 0) continue;

                // Interpolate edge vertices
                unsigned int edgeVerts[12] = {0};
                for (int e = 0; e < 12; ++e) {
                    if (!(mc::edgeTable[cubeIndex] & (1 << e))) continue;

                    auto key = std::make_tuple(ix * 3 + mc::edgeConnection[e][0],
                                                iy * 3 + mc::edgeConnection[e][0],
                                                iz, e);
                    // Unique key: cell position + edge index
                    key = std::make_tuple(ix, iy, iz, e);

                    auto it = edgeVertexMap.find(key);
                    if (it != edgeVertexMap.end()) {
                        edgeVerts[e] = it->second;
                    } else {
                        int c0 = mc::edgeConnection[e][0];
                        int c1 = mc::edgeConnection[e][1];
                        float v0 = cornerVals[c0];
                        float v1 = cornerVals[c1];
                        float t = (std::abs(v1 - v0) > 1e-10f) ? (isovalue - v0) / (v1 - v0) : 0.5f;
                        t = std::clamp(t, 0.0f, 1.0f);

                        float gx = ix + mc::cornerOffsets[c0][0] + t * (mc::cornerOffsets[c1][0] - mc::cornerOffsets[c0][0]);
                        float gy = iy + mc::cornerOffsets[c0][1] + t * (mc::cornerOffsets[c1][1] - mc::cornerOffsets[c0][1]);
                        float gz = iz + mc::cornerOffsets[c0][2] + t * (mc::cornerOffsets[c1][2] - mc::cornerOffsets[c0][2]);

                        // World position (mm -> world)
                        float wx = (minX + gx * voxelSize) * scale;
                        float wy = (minY + gy * voxelSize) * scale;
                        float wz = (minZ + gz * voxelSize) * scale;

                        // Normal from gradient
                        glm::vec3 n = computeGradient(gx, gy, gz);

                        unsigned int vi = (unsigned int)(vertices.size() / 6);
                        vertices.push_back(wx);
                        vertices.push_back(wy);
                        vertices.push_back(wz);
                        vertices.push_back(n.x);
                        vertices.push_back(n.y);
                        vertices.push_back(n.z);

                        edgeVertexMap[key] = vi;
                        edgeVerts[e] = vi;
                    }
                }

                // Generate triangles
                for (int t = 0; mc::triTable[cubeIndex][t] != -1; t += 3) {
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t]]);
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t + 1]]);
                    indices.push_back(edgeVerts[mc::triTable[cubeIndex][t + 2]]);
                }
            }
        }
    }

    if (vertices.empty() || indices.empty()) {
        // Even if marching cubes produces nothing, still add empty mesh
        StructureMesh mesh;
        auto color = structure->getColor();
        mesh.color = glm::vec3(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
        mesh.visible = structure->isVisible();
        mesh.indexCount = 0;
        mesh.center = glm::vec3(0.0f);
        mesh.vao = 0;
        mesh.vbo = 0;
        mesh.ebo = 0;
        m_meshes.push_back(mesh);
        Logger::info("    -> Skipped (marching cubes produced no geometry)");
        return;
    }

    // 7. Smooth normals by averaging face normals per vertex
    {
        size_t vertCount = vertices.size() / 6;
        std::vector<glm::vec3> smoothNormals(vertCount, glm::vec3(0.0f));

        for (size_t i = 0; i + 2 < indices.size(); i += 3) {
            unsigned int i0 = indices[i], i1 = indices[i+1], i2 = indices[i+2];
            glm::vec3 v0(vertices[i0*6], vertices[i0*6+1], vertices[i0*6+2]);
            glm::vec3 v1(vertices[i1*6], vertices[i1*6+1], vertices[i1*6+2]);
            glm::vec3 v2(vertices[i2*6], vertices[i2*6+1], vertices[i2*6+2]);
            glm::vec3 faceN = glm::cross(v1 - v0, v2 - v0);
            smoothNormals[i0] += faceN;
            smoothNormals[i1] += faceN;
            smoothNormals[i2] += faceN;
        }

        for (size_t i = 0; i < vertCount; ++i) {
            glm::vec3 n = glm::normalize(smoothNormals[i]);
            if (glm::length(smoothNormals[i]) < 1e-8f) n = glm::vec3(0,0,1);
            vertices[i*6+3] = n.x;
            vertices[i*6+4] = n.y;
            vertices[i*6+5] = n.z;
        }
    }

    // Compute bounding box center for sorting
    glm::vec3 minBound(1e30f), maxBound(-1e30f);
    size_t vertCount = vertices.size() / 6;
    for (size_t i = 0; i < vertCount; ++i) {
        glm::vec3 pos(vertices[i*6], vertices[i*6+1], vertices[i*6+2]);
        minBound = glm::min(minBound, pos);
        maxBound = glm::max(maxBound, pos);
    }
    glm::vec3 meshCenter = (minBound + maxBound) * 0.5f;

    // 8. Upload to GPU
    StructureMesh mesh;
    auto color = structure->getColor();
    mesh.color = glm::vec3(color[0] / 255.0f, color[1] / 255.0f, color[2] / 255.0f);
    mesh.visible = structure->isVisible();
    mesh.indexCount = indices.size();
    mesh.center = meshCenter;

    glGenVertexArrays(1, &mesh.vao);
    glGenBuffers(1, &mesh.vbo);
    glGenBuffers(1, &mesh.ebo);

    glBindVertexArray(mesh.vao);

    glBindBuffer(GL_ARRAY_BUFFER, mesh.vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh.ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(), GL_STATIC_DRAW);

    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    glBindVertexArray(0);

    m_meshes.push_back(mesh);
}

void StructureRenderer::buildMeshes() {
    std::lock_guard<std::mutex> lock(m_meshMutex);
    buildMeshes_unlocked();
}

void StructureRenderer::buildMeshes_unlocked() {
    // Called with m_meshMutex already held - no additional locking needed
    
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1, &mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1, &mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1, &mesh.ebo);
    }
    m_meshes.clear();

    if (!m_patientData || !m_patientData->getStructureSet()) return;

    auto* structures = m_patientData->getStructureSet();
    size_t structureCount = structures->getCount();
    
    Logger::info("Building meshes for " + std::to_string(structureCount) + " structures");
    
    for (size_t i = 0; i < structureCount; ++i) {
        const auto* structure = structures->getStructure(i);
        if (structure) {
            Logger::info("  Structure " + std::to_string(i) + ": " + structure->getName() + 
                        " (" + std::to_string(structure->getContourCount()) + " contours)");
        }
        tessellateStructure(i);
    }
    
    Logger::info("Built " + std::to_string(m_meshes.size()) + " meshes");
    
    // Verify alignment
    if (m_meshes.size() != structureCount) {
        Logger::error("MISMATCH: " + std::to_string(m_meshes.size()) + " meshes vs " + 
                     std::to_string(structureCount) + " structures!");
    }

    m_needsRebuild = false;
}

void StructureRenderer::render(const glm::mat4& view, const glm::mat4& projection) {
    std::lock_guard<std::mutex> lock(m_meshMutex);  // Synchronize with buildMeshes()
    
    if (!m_patientData || !m_patientData->getStructureSet()) return;

    if (m_needsRebuild) buildMeshes_unlocked();
    if (m_meshes.empty()) return;

    auto* structures = m_patientData->getStructureSet();
    size_t structureCount = structures->getCount();
    
    // Safety check: ensure mesh count matches structure count
    if (m_meshes.size() != structureCount) {
        Logger::warn("Mesh/structure count mismatch during render: " + 
                    std::to_string(m_meshes.size()) + " meshes vs " + 
                    std::to_string(structureCount) + " structures. Rebuilding...");
        buildMeshes();
        if (m_meshes.size() != structureCount) {
            Logger::error("Failed to fix mesh/structure alignment!");
            return;
        }
    }

    // -90° rotation to align with anatomical axes
    glm::mat4 model = glm::rotate(glm::mat4(1.0f), -glm::half_pi<float>(), glm::vec3(1.0f, 0.0f, 0.0f));

    // Extract camera position from inverse view matrix
    glm::mat4 invView = glm::inverse(view);
    glm::vec3 viewPos(invView[3][0], invView[3][1], invView[3][2]);

    // Two-light setup
    glm::vec3 lightDir1 = glm::normalize(glm::vec3(0.5f, 0.8f, 0.6f));
    glm::vec3 lightDir2 = glm::normalize(glm::vec3(-0.4f, -0.3f, 0.8f));

    // Enable blending and depth testing
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_TRUE); // Keep depth writes enabled
    glDisable(GL_CULL_FACE); // Two-sided rendering

    glUseProgram(m_shaderProgram);
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "model"), 1, GL_FALSE, glm::value_ptr(model));
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "view"), 1, GL_FALSE, glm::value_ptr(view));
    glUniformMatrix4fv(glGetUniformLocation(m_shaderProgram, "projection"), 1, GL_FALSE, glm::value_ptr(projection));
    glUniform3fv(glGetUniformLocation(m_shaderProgram, "viewPos"), 1, glm::value_ptr(viewPos));
    glUniform3fv(glGetUniformLocation(m_shaderProgram, "lightDir1"), 1, glm::value_ptr(lightDir1));
    glUniform3fv(glGetUniformLocation(m_shaderProgram, "lightDir2"), 1, glm::value_ptr(lightDir2));
    glUniform1f(glGetUniformLocation(m_shaderProgram, "opacity"), 0.30f);

    // Render all visible meshes
    for (size_t i = 0; i < m_meshes.size(); ++i) {
        // Skip empty meshes
        if (m_meshes[i].indexCount == 0) continue;
        
        // Check visibility
        const auto* structure = structures->getStructure(i);
        if (!structure || !structure->isVisible()) continue;
        
        const auto& mesh = m_meshes[i];
        glUniform3fv(glGetUniformLocation(m_shaderProgram, "structureColor"), 1, glm::value_ptr(mesh.color));

        glBindVertexArray(mesh.vao);
        glDrawElements(GL_TRIANGLES, (GLsizei)mesh.indexCount, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
    glUseProgram(0);
    glEnable(GL_CULL_FACE);
    glDisable(GL_BLEND);
}

void StructureRenderer::cleanup() {
    for (auto& mesh : m_meshes) {
        if (mesh.vao) glDeleteVertexArrays(1, &mesh.vao);
        if (mesh.vbo) glDeleteBuffers(1, &mesh.vbo);
        if (mesh.ebo) glDeleteBuffers(1, &mesh.ebo);
    }
    m_meshes.clear();

    if (m_shaderProgram) {
        glDeleteProgram(m_shaderProgram);
        m_shaderProgram = 0;
    }
}

} // namespace optirad
