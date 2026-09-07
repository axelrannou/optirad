#pragma once

#include "core/workflow/PlanAnalysis.hpp"
#include <fstream>
#include <string>
#include <vector>

namespace optirad {

/// Quote a CSV field if it contains a comma, quote, or newline.
inline std::string csvQuote(const std::string& field) {
    if (field.find_first_of(",\"\n") == std::string::npos) return field;
    std::string out = "\"";
    for (char c : field) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

/// Replace characters unsafe for filenames with '_'.
inline std::string sanitizeFilename(const std::string& name) {
    std::string out = name;
    for (char& c : out) {
        if (c == '/' || c == '\\' || c == ':' || c == '*' || c == '?' ||
            c == '"' || c == '<' || c == '>' || c == '|' || c == ' ')
            c = '_';
    }
    return out;
}

inline void writeDoseStatsCsvHeader(std::ofstream& f) {
    f << "Dose,Structure,Type,NumVoxels,Min(Gy),Max(Gy),Mean(Gy),Std(Gy),"
         "D2(Gy),D5(Gy),D50(Gy),D95(Gy),D98(Gy),V20(%),V40(%),V50(%),V60(%),CI,HI\n";
}

inline void writeDoseStatsCsvRows(std::ofstream& f, const std::string& doseName,
                                   const std::vector<StructureDoseStats>& stats) {
    for (const auto& s : stats) {
        f << csvQuote(doseName) << ',' << csvQuote(s.name) << ',' << csvQuote(s.type) << ','
          << s.numVoxels << ','
          << s.minDose << ',' << s.maxDose << ',' << s.meanDose << ',' << s.stdDose << ','
          << s.d2 << ',' << s.d5 << ',' << s.d50 << ',' << s.d95 << ',' << s.d98 << ','
          << s.v20 << ',' << s.v40 << ',' << s.v50 << ',' << s.v60 << ','
          << s.conformityIndex << ',' << s.homogeneityIndex << '\n';
    }
}

} // namespace optirad
