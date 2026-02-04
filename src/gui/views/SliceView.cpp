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
    m_patientData = data;
    
    if (m_patientData && m_patientData->getCTVolume()) {
        auto dims = m_patientData->getCTVolume()->getGrid().getDimensions();
        
        switch (m_orientation) {
            case SliceOrientation::Axial:
                m_maxSlice = dims[2] - 1;
                m_textureWidth = dims[0];
                m_textureHeight = dims[1];
                break;
            case SliceOrientation::Sagittal:
                m_maxSlice = dims[0] - 1;
                m_textureWidth = dims[1];
                m_textureHeight = dims[2];
                break;
            case SliceOrientation::Coronal:
                m_maxSlice = dims[1] - 1;
                m_textureWidth = dims[0];
                m_textureHeight = dims[2];
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
                    hu = ct->at(x, y, m_currentSlice);
                    break;
                case SliceOrientation::Sagittal:
                    hu = ct->at(m_currentSlice, x, y);
                    break;
                case SliceOrientation::Coronal:
                    hu = ct->at(x, m_currentSlice, y);
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
    
    // Calculate display size maintaining aspect ratio
    float aspect = static_cast<float>(m_textureWidth) / m_textureHeight;
    float displayWidth = avail.x;
    float displayHeight = avail.x / aspect;
    
    if (displayHeight > avail.y) {
        displayHeight = avail.y;
        displayWidth = avail.y * aspect;
    }
    
    // Center the image
    ImVec2 cursor = ImGui::GetCursorPos();
    ImGui::SetCursorPos(ImVec2(cursor.x + (avail.x - displayWidth) * 0.5f,
                               cursor.y + (avail.y - displayHeight) * 0.5f));
    
    // Display texture
    ImGui::Image(reinterpret_cast<void*>(static_cast<intptr_t>(m_textureID)),
                ImVec2(displayWidth, displayHeight));
}

void SliceView::update() {
    // No per-frame update needed
}

} // namespace optirad
