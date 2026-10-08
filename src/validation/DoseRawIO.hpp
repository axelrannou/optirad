#pragma once

#include "dose/DoseMatrix.hpp"
#include <memory>
#include <string>

namespace optirad {

/// Minimal binary dose format used to exchange dose cubes with other tools (e.g. matRad).
/// Layout: "ODOSE1\0\0", int32 n[3], double spacing[3], double origin[3] (LPS, mm), double data[n0*n1*n2].
/// Index order follows Grid: first index runs along patient y, second along x, third along z
/// (identical to a matRad [y,x,z] cube flattened column-major).
class DoseRawIO {
public:
    static void write(const std::string& path, const DoseMatrix& dose);
    static std::shared_ptr<DoseMatrix> read(const std::string& path);  // throws on error
};

} // namespace optirad
