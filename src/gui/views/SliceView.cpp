#include "SliceView.hpp"

namespace optirad {

void SliceView::update() {}

void SliceView::render() {
    // TODO: Render CT slice with dose overlay using OpenGL
    // 1. Extract slice from volume
    // 2. Upload to texture
    // 3. Render as textured quad
    // 4. Overlay dose with colormap
}

void SliceView::resize(int width, int height) {
    // TODO: Resize framebuffer
}

void SliceView::setCT(const CTVolume* ct) { m_ct = ct; }
void SliceView::setDose(const DoseVolume* dose) { m_dose = dose; }
void SliceView::setSliceIndex(int index) { m_sliceIndex = index; }
void SliceView::setOrientation(SliceOrientation orientation) { m_orientation = orientation; }

} // namespace optirad
