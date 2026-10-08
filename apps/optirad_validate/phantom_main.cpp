// Computes the water-phantom benchmark dose and writes it as an ODOSE1 file.
#include "validation/DoseRawIO.hpp"
#include "validation/WaterPhantom.hpp"

#include <algorithm>
#include <cstdlib>
#include <iostream>

using namespace optirad;

int main(int argc, char** argv) {
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: " << argv[0] << " <output.odose> [relative_threshold]\n";
        return 2;
    }

    try {
        WaterPhantomConfig cfg;
        if (argc == 3) cfg.relativeThreshold = std::atof(argv[2]);
        const auto r = WaterPhantom::compute(cfg);
        DoseRawIO::write(argv[1], *r.dose);

        // Depth is measured from the water surface.
        const auto pdd = sampleDepthDose(*r.dose, 0, 0, r.surfaceY, r.surfaceY + 250.0, 1.0);
        const auto it = std::max_element(pdd.dose.begin(), pdd.dose.end());
        const double dmax = *it;
        std::cout << "Water phantom, " << cfg.fieldSizeMm << " mm square field, SSD "
                  << 1000.0 - cfg.isocenterDepthMm << " mm\n"
                  << "  dmax depth: " << pdd.position[it - pdd.dose.begin()] - r.surfaceY << " mm\n";
        for (double depth : {50.0, 100.0, 200.0}) {
            const auto c = sampleDepthDose(*r.dose, 0, 0, r.surfaceY + depth, r.surfaceY + depth, 1.0);
            if (!c.dose.empty()) std::cout << "  PDD(" << depth << " mm): " << 100.0 * c.dose[0] / dmax << " %\n";
        }
        std::cout << "Written: " << argv[1] << "\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
