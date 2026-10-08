// Dose-grid and bixel-width convergence study on the water phantom.
// Study A: dose grid resolution vs a fine-grid reference (bixel width fixed).
// Study B: bixel width vs a fine-bixel reference (dose grid fixed).
// Each pairing is scored with 3%/2mm and 1%/1mm gamma for an open field and a stepped fluence.
#include "validation/GammaAnalysis.hpp"
#include "validation/WaterPhantom.hpp"

#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

using namespace optirad;

namespace {

constexpr double kHalf = 30.0;  // field half-size [mm]

WaterPhantomConfig baseConfig(double resolution, double bixelWidth) {
    WaterPhantomConfig c;
    c.resolutionMm = resolution;
    c.bixelWidthMm = bixelWidth;
    c.fieldSizeMm = 2 * kHalf;
    c.lateralHalfSizeMm = 90.0;
    c.depthMm = 210.0;
    c.isocenterDepthMm = 100.0;
    c.relativeThreshold = 1e-4;
    auto inside = [](double x, double z) { return std::abs(x) <= kHalf + 1e-6 && std::abs(z) <= kHalf + 1e-6; };
    c.fluences = {
        [=](double x, double z) { return inside(x, z) ? 1.0 : 0.0; },                       // open
        [=](double x, double z) { return inside(x, z) ? (x < 2.0 ? 1.0 : 0.4) : 0.0; }};    // step at x = 2 mm
    return c;
}

struct Row {
    std::string label;
    WaterPhantomResult result;
};

void score(const Row& ref, const Row& row, const char* fluenceName, size_t f) {
    GammaOptions g32; g32.doseDifferencePercent = 3.0; g32.distanceToAgreementMm = 2.0;
    GammaOptions g11; g11.doseDifferencePercent = 1.0; g11.distanceToAgreementMm = 1.0;
    const auto a = GammaAnalysis::compute(*ref.result.doses[f], *row.result.doses[f], g32);
    const auto b = GammaAnalysis::compute(*ref.result.doses[f], *row.result.doses[f], g11);
    std::printf("  %-14s %-6s  3%%/2mm: %6.2f%% (mean %.3f)   1%%/1mm: %6.2f%% (mean %.3f)   Dij %5.1fM nnz %6.1fs\n",
                row.label.c_str(), fluenceName, a.passRate, a.meanGamma, b.passRate, b.meanGamma,
                row.result.dijNonZeros / 1e6, row.result.dijSeconds);
}

void study(const char* title, const Row& ref, const std::vector<Row>& rows) {
    std::printf("\n%s (reference: %s, Dij %.1fM nnz, %.1fs)\n", title, ref.label.c_str(),
                ref.result.dijNonZeros / 1e6, ref.result.dijSeconds);
    for (const auto& row : rows) {
        score(ref, row, "open", 0);
        score(ref, row, "step", 1);
    }
}

} // namespace

int main() {
    try {
        {
            Row ref{"grid 1.25 mm", WaterPhantom::compute(baseConfig(1.25, 5.0))};
            std::vector<Row> rows;
            for (double res : {5.0, 3.0, 2.5, 2.0})
                rows.push_back({"grid " + std::to_string(res).substr(0, 4) + " mm", WaterPhantom::compute(baseConfig(res, 5.0))});
            study("Study A: dose grid resolution, bixel width 5 mm", ref, rows);
        }
        {
            Row ref{"bixel 1.5 mm", WaterPhantom::compute(baseConfig(2.5, 1.5))};
            std::vector<Row> rows;
            for (double bw : {7.0, 5.0, 3.0, 2.5})
                rows.push_back({"bixel " + std::to_string(bw).substr(0, 3) + " mm", WaterPhantom::compute(baseConfig(2.5, bw))});
            study("Study B: bixel width, dose grid 2.5 mm", ref, rows);
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
