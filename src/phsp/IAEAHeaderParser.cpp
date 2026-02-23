#include "phsp/IAEAHeaderParser.hpp"
#include "utils/Logger.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <cctype>

namespace optirad {

// ---------------------------------------------------------------------------
// Byte-order detection
// ---------------------------------------------------------------------------

bool IAEAHeaderInfo::needsByteSwap() const {
    // Detect host byte order
    uint32_t test = 1;
    bool hostLittleEndian = (*reinterpret_cast<uint8_t*>(&test) == 1);
    bool fileLittleEndian = (byteOrder == 1234);
    return hostLittleEndian != fileLittleEndian;
}

// ---------------------------------------------------------------------------
// Helper: trim whitespace
// ---------------------------------------------------------------------------

static std::string trim(const std::string& s) {
    size_t start = s.find_first_not_of(" \t\r\n");
    if (start == std::string::npos) return "";
    size_t end = s.find_last_not_of(" \t\r\n");
    return s.substr(start, end - start + 1);
}

// ---------------------------------------------------------------------------
// Helper: parse a single integer from a line (strips comments after //)
// ---------------------------------------------------------------------------

static int64_t parseIntValue(const std::string& line) {
    std::string clean = line;
    auto pos = clean.find("//");
    if (pos != std::string::npos) clean = clean.substr(0, pos);
    clean = trim(clean);
    if (clean.empty()) return 0;
    return std::stoll(clean);
}

static double parseDoubleValue(const std::string& line) {
    std::string clean = line;
    auto pos = clean.find("//");
    if (pos != std::string::npos) clean = clean.substr(0, pos);
    clean = trim(clean);
    if (clean.empty()) return 0.0;
    return std::stod(clean);
}

// ---------------------------------------------------------------------------
// IAEAHeaderParser::parse
// ---------------------------------------------------------------------------

IAEAHeaderInfo IAEAHeaderParser::parse(const std::string& headerFilePath) {
    Logger::debug("Parsing IAEA header: " + headerFilePath);

    std::ifstream file(headerFilePath);
    if (!file.is_open()) {
        throw std::runtime_error("IAEAHeaderParser: cannot open file: " + headerFilePath);
    }

    IAEAHeaderInfo info;
    std::string line;
    std::string currentKey;

    // Read all lines and dispatch based on $KEY: markers
    while (std::getline(file, line)) {
        std::string trimmed = trim(line);
        if (trimmed.empty()) continue;

        // Check for key marker
        if (trimmed[0] == '$' && trimmed.back() == ':') {
            currentKey = trimmed.substr(0, trimmed.size() - 1); // remove trailing ':'
            continue;
        }

        // If we were expecting the title, accumulate full line
        if (currentKey == "$IAEA_INDEX") {
            info.iaeaIndex = static_cast<int>(parseIntValue(trimmed));
            currentKey.clear();
        }
        else if (currentKey == "$TITLE") {
            info.title = trimmed;
            currentKey.clear();
        }
        else if (currentKey == "$FILE_TYPE") {
            info.fileType = static_cast<int>(parseIntValue(trimmed));
            currentKey.clear();
        }
        else if (currentKey == "$RECORD_CONTENTS") {
            // IAEA standard: 9 values with indices:
            //   [0]=Ek, [1]=X, [2]=Y, [3]=Z, [4]=U, [5]=V,
            //   [6]=Wt, [7]=ExtraFloats, [8]=ExtraLongs
            // Note: many header files have INCORRECT comments labeling these
            // as X,Y,Z,U,V,W but the actual IAEA spec maps index 0 to Energy.
            // W (direction cosine Z) is computed from U,V, never stored.
            static const int kNumFields = 9;
            int values[kNumFields];
            values[0] = static_cast<int>(parseIntValue(trimmed));
            for (int i = 1; i < kNumFields; ++i) {
                if (std::getline(file, line)) {
                    values[i] = static_cast<int>(parseIntValue(trim(line)));
                }
            }
            info.storesEnergy = (values[0] == 1);
            info.storesX = (values[1] == 1);
            info.storesY = (values[2] == 1);
            info.storesZ = (values[3] == 1);
            info.storesU = (values[4] == 1);
            info.storesV = (values[5] == 1);
            info.storesWeight = (values[6] == 1);
            info.numExtraFloats = values[7];
            info.numExtraLongs = values[8];
            currentKey.clear();
        }
        else if (currentKey == "$RECORD_CONSTANT") {
            info.constantWeight = parseDoubleValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$RECORD_LENGTH") {
            info.recordLength = static_cast<int>(parseIntValue(trimmed));
            currentKey.clear();
        }
        else if (currentKey == "$BYTE_ORDER") {
            info.byteOrder = static_cast<int>(parseIntValue(trimmed));
            currentKey.clear();
        }
        else if (currentKey == "$ORIG_HISTORIES") {
            info.origHistories = parseIntValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$PARTICLES") {
            info.totalParticles = parseIntValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$PHOTONS") {
            info.numPhotons = parseIntValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$ELECTRONS") {
            info.numElectrons = parseIntValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$POSITRONS") {
            info.numPositrons = parseIntValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$MONTE_CARLO_CODE_VERSION") {
            info.mcCodeVersion = trimmed;
            currentKey.clear();
        }
        else if (currentKey == "$GLOBAL_PHOTON_ENERGY_CUTOFF") {
            info.photonEnergyCutoff = parseDoubleValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$GLOBAL_PARTICLE_ENERGY_CUTOFF") {
            info.particleEnergyCutoff = parseDoubleValue(trimmed);
            currentKey.clear();
        }
        else if (currentKey == "$NOMINAL_SSD") {
            // Value may have "cm" suffix
            std::string val = trimmed;
            auto cmPos = val.find("cm");
            if (cmPos != std::string::npos) val = val.substr(0, cmPos);
            info.nominalSSD = parseDoubleValue(val);
            currentKey.clear();
        }
        else if (currentKey == "$COORDINATE_SYSTEM_DESCRIPTION") {
            // Multi-line description — accumulate until next $KEY
            info.coordinateSystemDesc = trimmed;
            // Peek ahead and accumulate non-key lines
            while (file.peek() != EOF) {
                auto pos = file.tellg();
                if (!std::getline(file, line)) break;
                std::string t = trim(line);
                if (!t.empty() && t[0] == '$') {
                    file.seekg(pos); // put line back
                    break;
                }
                info.coordinateSystemDesc += "\n" + t;
            }
            currentKey.clear();
        }
        else if (currentKey == "$STATISTICAL_INFORMATION_PARTICLES") {
            // Skip comment header line, then parse 3 particle lines
            // Format: weight wmin wmax <E> Emin Emax PARTICLE_TYPE
            for (int i = 0; i < 3; ++i) {
                if (!std::getline(file, line)) break;
                std::string t = trim(line);
                if (t.empty() || t[0] == '/') { --i; continue; } // skip comment lines

                std::istringstream iss(t);
                double weight, wmin, wmax, meanE, emin, emax;
                std::string ptype;
                if (iss >> weight >> wmin >> wmax >> meanE >> emin >> emax >> ptype) {
                    if (ptype == "PHOTONS") {
                        info.photonMeanEnergy = meanE;
                    } else if (ptype == "ELECTRONS") {
                        info.electronMeanEnergy = meanE;
                    } else if (ptype == "POSITRONS") {
                        info.positronMeanEnergy = meanE;
                    }
                }
            }
            currentKey.clear();
        }
        else if (currentKey == "$STATISTICAL_INFORMATION_GEOMETRY") {
            // 3 lines: x_range, y_range, z_range
            {
                std::istringstream iss(trimmed);
                iss >> info.xRange[0] >> info.xRange[1];
            }
            for (int i = 0; i < 2; ++i) {
                if (!std::getline(file, line)) break;
                std::string t = trim(line);
                if (t.empty()) { --i; continue; }
                std::istringstream iss(t);
                if (i == 0) iss >> info.yRange[0] >> info.yRange[1];
                else        iss >> info.zRange[0] >> info.zRange[1];
            }
            currentKey.clear();
        }
    }

    // Validate
    int expectedLen = info.computeRecordLength();
    if (expectedLen != info.recordLength) {
        Logger::warn("IAEAHeaderParser: computed record length (" + std::to_string(expectedLen) +
                     ") differs from declared (" + std::to_string(info.recordLength) + ")");
    }

    Logger::debug("IAEA header parsed: " + std::to_string(info.totalParticles) + " particles (" +
                  std::to_string(info.numPhotons) + " photons, " +
                  std::to_string(info.numElectrons) + " electrons, " +
                  std::to_string(info.numPositrons) + " positrons)");

    return info;
}

} // namespace optirad
