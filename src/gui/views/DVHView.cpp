#include "DVHView.hpp"

namespace optirad {

void DVHView::update() {}

void DVHView::render() {
    // TODO: Use ImPlot to render DVH curves
    // ImPlot::BeginPlot("DVH");
    // for (const auto& curve : m_curves) {
    //     ImPlot::SetNextLineStyle(ImVec4(curve.color[0], curve.color[1], curve.color[2], 1.0f));
    //     ImPlot::PlotLine(curve.structureName.c_str(), 
    //                      curve.doses.data(), curve.volumes.data(), curve.doses.size());
    // }
    // ImPlot::EndPlot();
}

void DVHView::resize(int width, int height) {}

void DVHView::addCurve(const DVHCurve& curve) { m_curves.push_back(curve); }
void DVHView::clearCurves() { m_curves.clear(); }

} // namespace optirad
