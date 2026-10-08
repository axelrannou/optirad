// Compare two dose files (DICOM RT Dose or .odose) with a 3D gamma analysis.
#include "io/DicomImporter.hpp"
#include "validation/DoseRawIO.hpp"
#include "validation/GammaAnalysis.hpp"

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

using namespace optirad;

static void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " <reference_rtdose.dcm> <evaluated_rtdose.dcm>"
                 " [--dd <percent>] [--dta <mm>] [--local] [--cutoff <percent>]\n";
}

static std::shared_ptr<DoseMatrix> loadDose(const std::string& path) {
    if (path.size() > 6 && path.compare(path.size() - 6, 6, ".odose") == 0) return DoseRawIO::read(path);
    DicomImporter importer;
    if (!importer.loadRTDose(path)) return nullptr;
    return importer.importRTDose().first;
}

int main(int argc, char** argv) {
    if (argc < 3) { usage(argv[0]); return 2; }

    GammaOptions opt;
    for (int i = 3; i < argc; ++i) {
        const bool hasValue = i + 1 < argc;
        if (!std::strcmp(argv[i], "--dd") && hasValue)          opt.doseDifferencePercent = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--dta") && hasValue)    opt.distanceToAgreementMm = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--cutoff") && hasValue) opt.lowDoseCutoffPercent = std::atof(argv[++i]);
        else if (!std::strcmp(argv[i], "--local"))              opt.globalNormalization = false;
        else { usage(argv[0]); return 2; }
    }

    std::shared_ptr<DoseMatrix> ref, eval;
    try {
        ref  = loadDose(argv[1]);
        eval = loadDose(argv[2]);
    } catch (const std::exception& e) {
        std::cerr << e.what() << "\n";
        return 1;
    }
    if (!ref || !eval) {
        std::cerr << "Failed to load dose from the given files\n";
        return 1;
    }

    GammaResult r;
    try {
        r = GammaAnalysis::compute(*ref, *eval, opt);
    } catch (const std::exception& e) {
        std::cerr << "Gamma analysis failed: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Gamma " << opt.doseDifferencePercent << "%/" << opt.distanceToAgreementMm << "mm ("
              << (opt.globalNormalization ? "global" : "local") << ", cutoff "
              << opt.lowDoseCutoffPercent << "%)\n"
              << "  Evaluated voxels: " << r.numEvaluated << "\n"
              << "  Pass rate:        " << r.passRate << " %\n"
              << "  Mean gamma:       " << r.meanGamma << "\n"
              << "  Gamma 95th pct:   " << r.gamma95 << "\n"
              << "  Max gamma:        " << r.maxGamma << "\n";
    return 0;
}
