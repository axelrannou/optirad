#pragma once

#include "IView.hpp"
#include "core/PatientData.hpp"
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
    void update() override;  // Declare but don't define inline
    void resize(int width, int height) override;
    std::string getName() const override;
    
    void setPatientData(PatientData* data);
    void setSliceIndex(size_t index);
    
private:
    void updateTexture();
    void renderControls();
    void renderSlice();
    
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
    bool m_needsUpdate = true;
};

} // namespace optirad
