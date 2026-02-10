#include "SliceView.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <GL/gl.h>
#include <algorithm>

namespace optirad {

SliceView::SliceView(SliceOrientation orientation) 
    : m_orientation(orientation) {
    glGenTextures(1, &m_textureID);
}

SliceView::~SliceView() {
    if (m_textureID) {
        glDeleteTextures(1, &m_textureID);
    }
}

std::string SliceView::getName() const {
    switch (m_orientation) {
        case SliceOrientation::Axial: return "Axial";
        case SliceOrientation::Sagittal: return "Sagittal";
        case SliceOrientation::Coronal: return "Coronal";
        default: return "Unknown";
    }
}

void SliceView::resize(int width, int height) {
    // Handle window resize if needed
}

void SliceView::setPatientData(PatientData* data) {
    // Only reset slices if this is new data
    if (m_patientData == data) {
        return;
    }
    
    m_patientData = data;
    
    if (m_patientData && m_patientData->getCTVolume()) {
        auto dims = m_patientData->getCTVolume()->getGrid().getDimensions();
        auto spacing = m_patientData->getCTVolume()->getGrid().getSpacing();
        
        switch (m_orientation) {
            case SliceOrientation::Axial:
                m_maxSlice = dims[2] - 1;
                m_textureWidth = dims[0];
                m_textureHeight = dims[1];
                m_physicalWidth = dims[0] * spacing[0];
                m_physicalHeight = dims[1] * spacing[1];
                break;
                
            case SliceOrientation::Sagittal:
                m_maxSlice = dims[0] - 1;
                m_textureWidth = dims[1];
                m_textureHeight = dims[2];
                m_physicalWidth = dims[1] * spacing[1];
                m_physicalHeight = dims[2] * spacing[2];
                break;
                
            case SliceOrientation::Coronal:
                m_maxSlice = dims[1] - 1;
                m_textureWidth = dims[0];
                m_textureHeight = dims[2];
                m_physicalWidth = dims[0] * spacing[0];
                m_physicalHeight = dims[2] * spacing[2];
                break;
        }
        
        m_currentSlice = m_maxSlice / 2;
        m_needsUpdate = true;
    }
}

void SliceView::setSliceIndex(size_t index) {
    if (index <= m_maxSlice) {
        m_currentSlice = index;
        m_needsUpdate = true;
    }
}

void SliceView::render() {
    std::string name = getName();
    ImGui::Begin(name.c_str());
    
    if (!m_patientData || !m_patientData->getCTVolume()) {
        ImGui::TextDisabled("No CT data loaded");
        ImGui::End();
        return;
    }
    
    renderControls();
    
    if (m_needsUpdate) {
        updateTexture();
        m_needsUpdate = false;
    }
    
    renderSlice();
    
    ImGui::End();
}

void SliceView::renderControls() {
    // Slice selector
    int sliceInt = static_cast<int>(m_currentSlice);
    if (ImGui::SliderInt("Slice", &sliceInt, 0, static_cast<int>(m_maxSlice))) {
        setSliceIndex(static_cast<size_t>(sliceInt));
    }
    
    // Show slice position in mm
    if (m_patientData && m_patientData->getCTVolume()) {
        const auto& grid = m_patientData->getCTVolume()->getGrid();
        auto origin = grid.getOrigin();
        auto spacing = grid.getSpacing();
        auto patientPos = grid.getPatientPosition();
        auto dims = grid.getDimensions();
        
        // Calculate position in patient coordinate system
        // DICOM uses LPS (Left, Posterior, Superior) coordinate system
        double slicePositionMM = 0.0;
        
        switch (m_orientation) {
            case SliceOrientation::Axial:
                // Z-axis: shows superior-inferior position
                // In LPS: larger Z = more superior (toward head)
                slicePositionMM = origin[2] + m_currentSlice * spacing[2];
                break;
                
            case SliceOrientation::Sagittal: {
                // X-axis: shows left-right position
                // Need to reverse: slice from left to right
                // Start from the last slice (most left) and go to right
                size_t reversedSlice = dims[0] - 1 - m_currentSlice;
                slicePositionMM = -(origin[0] + reversedSlice * spacing[0]);
                break;
            }
                
            case SliceOrientation::Coronal: {
                // Y-axis: shows anterior-posterior position
                // Need to reverse: slice from back to front (posterior to anterior)
                size_t reversedSlice = dims[1] - 1 - m_currentSlice;
                slicePositionMM = -(origin[1] + reversedSlice * spacing[1]);
                break;
            }
        }
        
        ImGui::SameLine();
        ImGui::Text("(%.1f mm)", slicePositionMM);
        
        // Debug tooltip showing both LPS and RAS coordinates
        if (ImGui::IsItemHovered()) {
            ImGui::BeginTooltip();
            
            switch (m_orientation) {
                case SliceOrientation::Axial: {
                    double lpsZ = origin[2] + m_currentSlice * spacing[2];
                    ImGui::Text("Patient Position: %s", patientPos.c_str());
                    ImGui::Separator();
                    ImGui::Text("Z (LPS Superior): %.1f mm", lpsZ);
                    ImGui::Text("Slice %zu/%zu", m_currentSlice, m_maxSlice);
                    break;
                }
                case SliceOrientation::Sagittal: {
                    size_t reversedSlice = dims[0] - 1 - m_currentSlice;
                    double lpsX = origin[0] + reversedSlice * spacing[0];
                    ImGui::Text("Patient Position: %s", patientPos.c_str());
                    ImGui::Separator();
                    ImGui::Text("X (LPS Left): %.1f mm", lpsX);
                    ImGui::Text("R (RAS Right): %.1f mm", -lpsX);
                    ImGui::Text("Slice %zu/%zu (reversed: %zu)", m_currentSlice, m_maxSlice, reversedSlice);
                    break;
                }
                case SliceOrientation::Coronal: {
                    size_t reversedSlice = dims[1] - 1 - m_currentSlice;
                    double lpsY = origin[1] + reversedSlice * spacing[1];
                    ImGui::Text("Patient Position: %s", patientPos.c_str());
                    ImGui::Separator();
                    ImGui::Text("Y (LPS Posterior): %.1f mm", lpsY);
                    ImGui::Text("A (RAS Anterior): %.1f mm", -lpsY);
                    ImGui::Text("Slice %zu/%zu (reversed: %zu)", m_currentSlice, m_maxSlice, reversedSlice);
                    break;
                }
            }
            
            ImGui::EndTooltip();
        }
    }
    
    // Contour display toggle
    ImGui::Checkbox("Show Contours", &m_showContours);
    ImGui::SameLine();
    if (m_showContours) {
        ImGui::SetNextItemWidth(100);
        ImGui::SliderFloat("Contour Thickness", &m_contourThickness, 0.5f, 5.0f);
    }
    
    // Window/Level controls
    if (ImGui::SliderInt("Window Width", &m_windowWidth, 1, 2000)) {
        m_needsUpdate = true;
    }
    if (ImGui::SliderInt("Window Center", &m_windowCenter, -1000, 1000)) {
        m_needsUpdate = true;
    }
    
    // Preset buttons
    if (ImGui::Button("Soft Tissue")) {
        m_windowCenter = 40;
        m_windowWidth = 400;
        m_needsUpdate = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Lung")) {
        m_windowCenter = -600;
        m_windowWidth = 1500;
        m_needsUpdate = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Bone")) {
        m_windowCenter = 400;
        m_windowWidth = 1800;
        m_needsUpdate = true;
    }
    
    // Show geometric info tooltip
    if (ImGui::IsItemHovered() || ImGui::IsWindowHovered()) {
        if (m_patientData && m_patientData->getCTVolume()) {
            const auto& grid = m_patientData->getCTVolume()->getGrid();
            
            if (ImGui::BeginTooltip()) {
                ImGui::Text("Patient Position: %s", grid.getPatientPosition().c_str());
                ImGui::Text("Slice Thickness: %.2f mm", grid.getSliceThickness());
                
                auto spacing = grid.getSpacing();
                ImGui::Text("Actual Spacing: %.2f x %.2f x %.2f mm", 
                           spacing[0], spacing[1], spacing[2]);
                
                ImGui::EndTooltip();
            }
        }
    }
    
    ImGui::Separator();
}

void SliceView::updateTexture() {
    auto* ct = m_patientData->getCTVolume();
    if (!ct) return;
    
    auto dims = ct->getGrid().getDimensions();
    
    // Extract slice data and apply window/level
    std::vector<unsigned char> pixels(m_textureWidth * m_textureHeight * 3);
    
    int minHU = m_windowCenter - m_windowWidth / 2;
    int maxHU = m_windowCenter + m_windowWidth / 2;
    
    for (int y = 0; y < m_textureHeight; ++y) {
        for (int x = 0; x < m_textureWidth; ++x) {
            int16_t hu = 0;
            
            // Get HU value based on orientation
            switch (m_orientation) {
                case SliceOrientation::Axial:
                    // Axial: Standard top-down view
                    hu = ct->at(x, y, m_currentSlice);
                    break;
                    
                case SliceOrientation::Sagittal:
                    // Sagittal: Y-Z plane, flip Y (Z increases superior, flip for screen)
                    {
                        int flippedY = m_textureHeight - 1 - y;
                        hu = ct->at(m_currentSlice, x, flippedY);
                    }
                    break;
                    
                case SliceOrientation::Coronal:
                    // Coronal: X-Z plane, flip both Y (for superior direction) and slice direction (posterior to anterior)
                    {
                        int flippedY = m_textureHeight - 1 - y;
                        size_t reversedSlice = dims[1] - 1 - m_currentSlice;
                        hu = ct->at(x, reversedSlice, flippedY);
                    }
                    break;
            }
            
            // Apply window/level
            unsigned char intensity = 0;
            if (hu <= minHU) {
                intensity = 0;
            } else if (hu >= maxHU) {
                intensity = 255;
            } else {
                intensity = static_cast<unsigned char>(255.0 * (hu - minHU) / (maxHU - minHU));
            }
            
            // Grayscale RGB
            size_t idx = 3 * (y * m_textureWidth + x);
            pixels[idx + 0] = intensity;
            pixels[idx + 1] = intensity;
            pixels[idx + 2] = intensity;
        }
    }
    
    // Upload to OpenGL texture
    glBindTexture(GL_TEXTURE_2D, m_textureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_textureWidth, m_textureHeight, 
                 0, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void SliceView::renderSlice() {
    // Get available space
    ImVec2 avail = ImGui::GetContentRegionAvail();
    
    // Calculate display size using PHYSICAL aspect ratio (like 3D Slicer)
    float physicalAspect = static_cast<float>(m_physicalWidth / m_physicalHeight);
    
    float displayWidth = avail.x;
    float displayHeight = avail.x / physicalAspect;
    
    if (displayHeight > avail.y) {
        displayHeight = avail.y;
        displayWidth = avail.y * physicalAspect;
    }
    
    // Center the image
    ImVec2 cursor = ImGui::GetCursorPos();
    ImVec2 imagePos = ImVec2(cursor.x + (avail.x - displayWidth) * 0.5f,
                             cursor.y + (avail.y - displayHeight) * 0.5f);
    ImGui::SetCursorPos(imagePos);
    
    // Store image bounds for contour rendering
    ImVec2 imageMin = ImGui::GetCursorScreenPos();
    ImVec2 imageMax = ImVec2(imageMin.x + displayWidth, imageMin.y + displayHeight);
    
    // Display texture
    ImGui::Image(reinterpret_cast<void*>(static_cast<intptr_t>(m_textureID)),
                ImVec2(displayWidth, displayHeight));
    
    // Overlay contours if enabled
    if (m_showContours && m_patientData && m_patientData->getStructureSet()) {
        renderContours();
    }
}

void SliceView::renderContours()
{
    // Only render contours for Axial view
    // RT-STRUCT contours are natively axial (constant Z per contour)
    if (m_orientation != SliceOrientation::Axial) {
        return;
    }
    
    auto* structures = m_patientData->getStructureSet();
    if (!structures) return;

    const auto& grid = m_patientData->getCTVolume()->getGrid();
    auto dims = grid.getDimensions();

    ImVec2 imageMin = ImGui::GetItemRectMin();
    ImVec2 imageMax = ImGui::GetItemRectMax();
    ImVec2 imageSize = ImVec2(imageMax.x - imageMin.x, imageMax.y - imageMin.y);

    ImDrawList* drawList = ImGui::GetWindowDrawList();

    double sliceIndexVoxel = (double)m_currentSlice;
    constexpr double sliceTol = 0.5;

    for (size_t si = 0; si < structures->getCount(); ++si) {
        const auto* structure = structures->getStructure(si);
        if (!structure || !structure->isVisible()) continue;

        auto c = structure->getColor();
        ImU32 imColor = IM_COL32(c[0], c[1], c[2], 255);

        for (const auto& contour : structure->getContours()) {
            if (contour.points.size() < 2) continue;

            // Convert contour points to voxel space
            std::vector<Vec3> vox;
            vox.reserve(contour.points.size());
            for (const auto& p : contour.points) {
                vox.push_back(grid.patientToVoxel({p[0], p[1], p[2]}));
            }

            // Check if contour is on current slice
            double z = vox[0][2];
            if (std::abs(z - sliceIndexVoxel) > sliceTol)
                continue;

            // Project to screen coordinates
            std::vector<ImVec2> poly;
            for (const auto& v : vox) {
                double u = v[0] / dims[0];
                double w = v[1] / dims[1];

                poly.push_back({
                    imageMin.x + (float)(u * imageSize.x),
                    imageMin.y + (float)(w * imageSize.y)
                });
            }

            if (poly.size() >= 2) {
                drawList->AddPolyline(
                    poly.data(),
                    (int)poly.size(),
                    imColor,
                    ImDrawFlags_Closed,
                    m_contourThickness
                );
            }
        }
    }
}

void SliceView::update() {
    // No per-frame update needed
}

} // namespace optirad
