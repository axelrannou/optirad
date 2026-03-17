#pragma once

#include "IView.hpp"
#include "core/PatientData.hpp"
#include "dose/DoseMatrix.hpp"
#include "geometry/Grid.hpp"
#include <memory>
#include <string>

namespace optirad {

enum class SliceOrientation {
    Axial,      // X-Y plane (Z slices)
    Sagittal,   // Y-Z plane (X slices)
    Coronal     // X-Z plane (Y slices)
};

class SliceView : public IView {
public:
    SliceView(SliceOrientation orientation);
    ~SliceView();
    
    void render() override;
    void update() override;
    void resize(int width, int height) override;
    std::string getName() const override;
    
    void setPatientData(PatientData* data);
    void setSliceIndex(size_t index);

    /// Set dose overlay data (nullptr to disable)
    void setDoseData(const DoseMatrix* dose, const Grid* doseGrid);
    
private:
    void updateTexture();
    void updateDoseTexture();
    void renderControls();
    void renderSlice();
    void renderContours();
    
    /// Jet colormap: value in [0,1] → R,G,B
    static void jetColormap(float t, unsigned char& r, unsigned char& g, unsigned char& b);
    
    PatientData* m_patientData = nullptr;
    SliceOrientation m_orientation;
    
    size_t m_currentSlice = 0;
    size_t m_maxSlice = 0;
    
    // Window/Level for CT display
    int m_windowCenter = 40;   // HU
    int m_windowWidth = 400;   // HU
    
    // OpenGL texture
    unsigned int m_textureID = 0;
    int m_textureWidth = 0;
    int m_textureHeight = 0;
    int m_lastTextureWidth = -1;
    int m_lastTextureHeight = -1;
    bool m_needsUpdate = true;
    
    // Physical dimensions for aspect ratio correction
    double m_physicalWidth = 1.0;   // mm
    double m_physicalHeight = 1.0;  // mm
    
    // Contour display settings
    float m_contourThickness = 2.0f;
    bool m_showContours = true;

    // Zoom and pan (UV-based)
    float m_zoom = 1.0f;
    float m_panU = 0.5f;   // UV center X [0,1]
    float m_panV = 0.5f;   // UV center Y [0,1]

    // Rendering state (set during renderSlice, used by renderContours)
    float m_visUvMinX = 0.0f, m_visUvMinY = 0.0f;
    float m_visUvMaxX = 1.0f, m_visUvMaxY = 1.0f;
    float m_imgScreenMinX = 0.0f, m_imgScreenMinY = 0.0f;
    float m_imgScreenMaxX = 0.0f, m_imgScreenMaxY = 0.0f;

    // Dose overlay
    const DoseMatrix* m_doseData = nullptr;
    const Grid* m_doseGrid = nullptr;
    unsigned int m_doseTextureID = 0;
    bool m_showDose = true;
    float m_doseAlpha = 0.8f;
    float m_doseThresholdPct = 0.5f; // % of max dose below which dose is not shown
    bool m_doseNeedsUpdate = true;
};

} // namespace optirad
