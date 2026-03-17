#include "SliceView.hpp"
#include "utils/Logger.hpp"
#include <imgui.h>
#include <GL/gl.h>
#include <algorithm>

namespace optirad {

SliceView::SliceView(SliceOrientation orientation) 
    : m_orientation(orientation) {
    glGenTextures(1, &m_textureID);
    glGenTextures(1, &m_doseTextureID);
}

SliceView::~SliceView() {
    if (m_textureID) {
        glDeleteTextures(1, &m_textureID);
    }
    if (m_doseTextureID) {
        glDeleteTextures(1, &m_doseTextureID);
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
        m_doseNeedsUpdate = true;
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
    // Line 1: Compact slice selector + position
    int sliceInt = static_cast<int>(m_currentSlice);
    ImGui::SetNextItemWidth(ImGui::GetContentRegionAvail().x - 70);
    if (ImGui::SliderInt("##Slice", &sliceInt, 0, static_cast<int>(m_maxSlice))) {
        setSliceIndex(static_cast<size_t>(sliceInt));
    }
    
    if (m_patientData && m_patientData->getCTVolume()) {
        const auto& grid = m_patientData->getCTVolume()->getGrid();
        auto origin = grid.getOrigin();
        auto spacing = grid.getSpacing();
        auto dims = grid.getDimensions();
        
        double slicePositionMM = 0.0;
        switch (m_orientation) {
            case SliceOrientation::Axial:
                slicePositionMM = origin[2] + m_currentSlice * spacing[2];
                break;
            case SliceOrientation::Sagittal: {
                size_t reversedSlice = dims[0] - 1 - m_currentSlice;
                slicePositionMM = -(origin[0] + reversedSlice * spacing[0]);
                break;
            }
            case SliceOrientation::Coronal: {
                size_t reversedSlice = dims[1] - 1 - m_currentSlice;
                slicePositionMM = -(origin[1] + reversedSlice * spacing[1]);
                break;
            }
        }
        ImGui::SameLine();
        ImGui::Text("%.1f", slicePositionMM);
    }
    
    // Line 2: Show Contours + Show Dose on same line
    ImGui::Checkbox("Contours", &m_showContours);
    if (m_doseData) {
        ImGui::SameLine();
        ImGui::Checkbox("Dose", &m_showDose);
    }
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
    
    // Only reallocate texture if dimensions changed
    if (m_textureWidth != m_lastTextureWidth || m_textureHeight != m_lastTextureHeight) {
        // Allocate new texture memory
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, m_textureWidth, m_textureHeight, 
                     0, GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        
        m_lastTextureWidth = m_textureWidth;
        m_lastTextureHeight = m_textureHeight;
    } else {
        // Update existing texture (more efficient)
        glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, m_textureWidth, m_textureHeight,
                        GL_RGB, GL_UNSIGNED_BYTE, pixels.data());
    }
    
    glBindTexture(GL_TEXTURE_2D, 0);
}

void SliceView::renderSlice() {
    ImVec2 avail = ImGui::GetContentRegionAvail();
    
    // Reserve space for preset buttons on the right
    float buttonW = 28.0f;
    float gap = 4.0f;
    float imgAreaW = avail.x - buttonW - gap;
    float imgAreaH = avail.y;
    
    if (imgAreaW < 10 || imgAreaH < 10) return;
    
    // Calculate display size using PHYSICAL aspect ratio
    float physicalAspect = static_cast<float>(m_physicalWidth / m_physicalHeight);
    
    float displayWidth = imgAreaW;
    float displayHeight = imgAreaW / physicalAspect;
    
    if (displayHeight > imgAreaH) {
        displayHeight = imgAreaH;
        displayWidth = imgAreaH * physicalAspect;
    }
    
    // Center the image in the image area
    ImVec2 cursor = ImGui::GetCursorPos();
    ImVec2 imagePos = ImVec2(cursor.x + (imgAreaW - displayWidth) * 0.5f,
                             cursor.y + (imgAreaH - displayHeight) * 0.5f);
    
    // Calculate UV coordinates based on zoom and pan
    float halfU = 0.5f / m_zoom;
    float halfV = 0.5f / m_zoom;
    ImVec2 uv0(m_panU - halfU, m_panV - halfV);
    ImVec2 uv1(m_panU + halfU, m_panV + halfV);
    
    // Store for contour rendering
    m_visUvMinX = uv0.x; m_visUvMinY = uv0.y;
    m_visUvMaxX = uv1.x; m_visUvMaxY = uv1.y;
    
    ImGui::SetCursorPos(imagePos);
    
    // Store image screen bounds for contour rendering
    ImVec2 imageScreenMin = ImGui::GetCursorScreenPos();
    m_imgScreenMinX = imageScreenMin.x;
    m_imgScreenMinY = imageScreenMin.y;
    m_imgScreenMaxX = imageScreenMin.x + displayWidth;
    m_imgScreenMaxY = imageScreenMin.y + displayHeight;
    
    // Display CT texture with zoom UVs
    ImGui::Image(reinterpret_cast<void*>(static_cast<intptr_t>(m_textureID)),
                ImVec2(displayWidth, displayHeight), uv0, uv1);
    
    // Check hover state on CT image
    bool imageHovered = ImGui::IsItemHovered();
    
    // Dose overlay with same UVs
    if (m_showDose && m_doseData && m_doseGrid) {
        if (m_doseNeedsUpdate) {
            updateDoseTexture();
            m_doseNeedsUpdate = false;
        }
        ImGui::SetCursorPos(imagePos);
        ImGui::Image(reinterpret_cast<void*>(static_cast<intptr_t>(m_doseTextureID)),
                    ImVec2(displayWidth, displayHeight), uv0, uv1);
    }
    
    // Overlay contours if enabled
    if (m_showContours && m_patientData && m_patientData->getStructureSet()) {
        renderContours();
    }
    
    // --- Mouse interactions ---
    
    // Left-click + drag vertical = zoom (like 3D Slicer)
    if (imageHovered && ImGui::IsMouseDragging(ImGuiMouseButton_Left)) {
        float deltaY = ImGui::GetIO().MouseDelta.y;
        m_zoom *= 1.0f - deltaY * 0.005f;
        m_zoom = std::clamp(m_zoom, 0.5f, 20.0f);
    }
    
    // Middle-click + drag = pan (when zoomed)
    if (imageHovered && ImGui::IsMouseDragging(ImGuiMouseButton_Middle)) {
        ImVec2 delta = ImGui::GetIO().MouseDelta;
        float uvPerPixelX = (uv1.x - uv0.x) / displayWidth;
        float uvPerPixelY = (uv1.y - uv0.y) / displayHeight;
        m_panU -= delta.x * uvPerPixelX;
        m_panV -= delta.y * uvPerPixelY;
        m_panU = std::clamp(m_panU, 0.0f, 1.0f);
        m_panV = std::clamp(m_panV, 0.0f, 1.0f);
    }
    
    // Scroll wheel = change slice
    if (imageHovered) {
        float wheel = ImGui::GetIO().MouseWheel;
        if (wheel != 0.0f) {
            int newSlice = static_cast<int>(m_currentSlice) - static_cast<int>(wheel);
            newSlice = std::clamp(newSlice, 0, static_cast<int>(m_maxSlice));
            setSliceIndex(static_cast<size_t>(newSlice));
        }
    }
    
    // Right-click context menu for W/L and settings
    if (imageHovered && ImGui::IsMouseClicked(ImGuiMouseButton_Right)) {
        ImGui::OpenPopup("##SliceCtx");
    }
    if (ImGui::BeginPopup("##SliceCtx")) {
        ImGui::Text("Window / Level");
        ImGui::Separator();
        ImGui::SetNextItemWidth(150);
        if (ImGui::SliderInt("Width", &m_windowWidth, 1, 2000)) m_needsUpdate = true;
        ImGui::SetNextItemWidth(150);
        if (ImGui::SliderInt("Center", &m_windowCenter, -1000, 1000)) m_needsUpdate = true;
        ImGui::Separator();
        if (ImGui::MenuItem("Soft Tissue (W:400 C:40)"))    { m_windowCenter = 40; m_windowWidth = 400; m_needsUpdate = true; }
        if (ImGui::MenuItem("Lung (W:1500 C:-600)"))        { m_windowCenter = -600; m_windowWidth = 1500; m_needsUpdate = true; }
        if (ImGui::MenuItem("Bone (W:1800 C:400)"))         { m_windowCenter = 400; m_windowWidth = 1800; m_needsUpdate = true; }
        if (m_showContours) {
            ImGui::Separator();
            ImGui::SetNextItemWidth(120);
            ImGui::SliderFloat("Contour Size", &m_contourThickness, 0.5f, 5.0f);
        }
        if (m_showDose && m_doseData) {
            ImGui::Separator();
            ImGui::SetNextItemWidth(120);
            if (ImGui::SliderFloat("Dose Opacity", &m_doseAlpha, 0.0f, 1.0f)) m_doseNeedsUpdate = true;
            ImGui::SetNextItemWidth(120);
            if (ImGui::SliderFloat("Dose Thresh%", &m_doseThresholdPct, 0.0f, 50.0f)) m_doseNeedsUpdate = true;
        }
        ImGui::Separator();
        if (ImGui::MenuItem("Reset Zoom")) { m_zoom = 1.0f; m_panU = 0.5f; m_panV = 0.5f; }
        ImGui::EndPopup();
    }
    
    // --- Preset buttons vertically on the right ---
    float btnH = ImGui::GetFrameHeight();
    float btnSpacing = ImGui::GetStyle().ItemSpacing.y;
    float totalBtnH = btnH * 3 + btnSpacing * 2;
    float btnStartY = cursor.y + (imgAreaH - totalBtnH) * 0.5f;
    float btnX = cursor.x + imgAreaW + gap;
    
    ImGui::SetCursorPos(ImVec2(btnX, btnStartY));
    if (ImGui::Button("ST", ImVec2(buttonW, 0))) {
        m_windowCenter = 40; m_windowWidth = 400; m_needsUpdate = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Soft Tissue (W:400 C:40)");
    ImGui::SetCursorPosX(btnX);
    if (ImGui::Button("Lg", ImVec2(buttonW, 0))) {
        m_windowCenter = -600; m_windowWidth = 1500; m_needsUpdate = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Lung (W:1500 C:-600)");
    ImGui::SetCursorPosX(btnX);
    if (ImGui::Button("Bn", ImVec2(buttonW, 0))) {
        m_windowCenter = 400; m_windowWidth = 1800; m_needsUpdate = true;
    }
    if (ImGui::IsItemHovered()) ImGui::SetTooltip("Bone (W:1800 C:400)");
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

    ImVec2 imageMin(m_imgScreenMinX, m_imgScreenMinY);
    ImVec2 imageMax(m_imgScreenMaxX, m_imgScreenMaxY);
    ImVec2 imageSize = ImVec2(imageMax.x - imageMin.x, imageMax.y - imageMin.y);

    ImDrawList* drawList = ImGui::GetWindowDrawList();
    
    // Clip contours to image bounds
    drawList->PushClipRect(imageMin, imageMax, true);

    double sliceIndexVoxel = (double)m_currentSlice;
    constexpr double sliceTol = 0.5;
    
    float uvRangeX = m_visUvMaxX - m_visUvMinX;
    float uvRangeY = m_visUvMaxY - m_visUvMinY;
    if (uvRangeX <= 0.0f || uvRangeY <= 0.0f) { drawList->PopClipRect(); return; }

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

            // Project to screen coordinates (zoom-aware via UV mapping)
            std::vector<ImVec2> poly;
            for (const auto& v : vox) {
                double u = v[1] / dims[1];
                double w = v[0] / dims[0];

                // Map from full UV [0,1] to visible UV range with zoom
                float screenX = imageMin.x + (float)((u - m_visUvMinX) / uvRangeX) * imageSize.x;
                float screenY = imageMin.y + (float)((w - m_visUvMinY) / uvRangeY) * imageSize.y;

                poly.push_back({screenX, screenY});
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
    
    drawList->PopClipRect();
}

void SliceView::update() {
    // No per-frame update needed
}

void SliceView::setDoseData(const DoseMatrix* dose, const Grid* doseGrid) {
    const bool doseSourceChanged = (m_doseData != dose) || (m_doseGrid != doseGrid);
    m_doseData = dose;
    m_doseGrid = doseGrid;
    if (doseSourceChanged) {
        m_doseNeedsUpdate = true;
    }
}

void SliceView::jetColormap(float t, unsigned char& r, unsigned char& g, unsigned char& b) {
    // Standard jet colormap: blue → cyan → green → yellow → red
    t = std::clamp(t, 0.0f, 1.0f);
    float r4 = std::clamp(1.5f - std::abs(t - 0.75f) * 4.0f, 0.0f, 1.0f);
    float g4 = std::clamp(1.5f - std::abs(t - 0.50f) * 4.0f, 0.0f, 1.0f);
    float b4 = std::clamp(1.5f - std::abs(t - 0.25f) * 4.0f, 0.0f, 1.0f);
    r = static_cast<unsigned char>(r4 * 255);
    g = static_cast<unsigned char>(g4 * 255);
    b = static_cast<unsigned char>(b4 * 255);
}

void SliceView::updateDoseTexture() {
    if (!m_doseData || !m_doseGrid || !m_patientData || !m_patientData->getCTVolume()) return;

    const auto& ctGrid = m_patientData->getCTVolume()->getGrid();
    auto ctDims = ctGrid.getDimensions();
    auto ctSpacing = ctGrid.getSpacing();
    auto ctOrigin = ctGrid.getOrigin();

    auto doseDims = m_doseGrid->getDimensions();
    auto doseSpacing = m_doseGrid->getSpacing();
    auto doseOrigin = m_doseGrid->getOrigin();

    double maxDose = m_doseData->getMax();
    if (maxDose <= 0.0) return;
    double threshold = maxDose * (m_doseThresholdPct / 100.0f);

    // Create RGBA pixels on CT grid dimensions (same as CT texture)
    std::vector<unsigned char> pixels(m_textureWidth * m_textureHeight * 4, 0);

    for (int y = 0; y < m_textureHeight; ++y) {
        for (int x = 0; x < m_textureWidth; ++x) {
            // Map CT texture pixel (x,y) + current slice to CT voxel coordinates.
            // Note: Volume::at(i,j,k) with dims[0]=rows stores data as
            // data[i + j*rows], so at(x,y,slice) for a square image gives
            // DICOM pixel(row=y, col=x). To get the correct LPS position via
            // voxelToPatient (where i→colDir=row direction, j→rowDir=col direction),
            // we must swap i↔j to match the actual DICOM row/column being displayed.
            double ctI = 0, ctJ = 0, ctK = 0;
            switch (m_orientation) {
                case SliceOrientation::Axial:
                    ctI = y; ctJ = x; ctK = m_currentSlice;
                    break;
                case SliceOrientation::Sagittal:
                    ctI = x; ctJ = m_currentSlice; ctK = m_textureHeight - 1 - y;
                    break;
                case SliceOrientation::Coronal:
                    ctI = ctDims[1] - 1 - m_currentSlice; ctJ = x; ctK = m_textureHeight - 1 - y;
                    break;
            }

            // Convert CT voxel to patient LPS
            Vec3 lps = ctGrid.voxelToPatient({ctI, ctJ, ctK});

            // Convert LPS to dose grid voxel
            Vec3 doseVox = m_doseGrid->patientToVoxel(lps);
            int di = static_cast<int>(std::round(doseVox[0]));
            int dj = static_cast<int>(std::round(doseVox[1]));
            int dk = static_cast<int>(std::round(doseVox[2]));

            if (di < 0 || di >= (int)doseDims[0] ||
                dj < 0 || dj >= (int)doseDims[1] ||
                dk < 0 || dk >= (int)doseDims[2]) continue;

            double dose = m_doseData->at(di, dj, dk);
            if (dose < threshold) continue;

            float t = static_cast<float>(dose / maxDose);
            unsigned char r, g, b;
            jetColormap(t, r, g, b);

            size_t idx = 4 * (y * m_textureWidth + x);
            pixels[idx + 0] = r;
            pixels[idx + 1] = g;
            pixels[idx + 2] = b;
            pixels[idx + 3] = static_cast<unsigned char>(m_doseAlpha * 255);
        }
    }

    glBindTexture(GL_TEXTURE_2D, m_doseTextureID);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, m_textureWidth, m_textureHeight,
                 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glBindTexture(GL_TEXTURE_2D, 0);
}

} // namespace optirad
