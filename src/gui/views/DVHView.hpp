#pragma once

#include "IView.hpp"
#include <vector>

namespace optirad {

struct DVHCurve {
    std::string structureName;
    std::vector<double> doses;      // x-axis
    std::vector<double> volumes;    // y-axis (0-100%)
    float color[3];
};

class DVHView : public IView {
public:
    std::string getName() const override { return "DVH"; }
    void update() override;
    void render() override;
    void resize(int width, int height) override;

    void addCurve(const DVHCurve& curve);
    void clearCurves();

private:
    std::vector<DVHCurve> m_curves;
};

} // namespace optirad
