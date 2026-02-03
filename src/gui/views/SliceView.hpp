#pragma once

#include "IView.hpp"
#include "Volume.hpp"

namespace optirad {

enum class SliceOrientation { Axial, Sagittal, Coronal };

class SliceView : public IView {
public:
    std::string getName() const override { return "Slice View"; }
    void update() override;
    void render() override;
    void resize(int width, int height) override;

    void setCT(const CTVolume* ct);
    void setDose(const DoseVolume* dose);
    void setSliceIndex(int index);
    void setOrientation(SliceOrientation orientation);

private:
    const CTVolume* m_ct = nullptr;
    const DoseVolume* m_dose = nullptr;
    int m_sliceIndex = 0;
    SliceOrientation m_orientation = SliceOrientation::Axial;
};

} // namespace optirad
