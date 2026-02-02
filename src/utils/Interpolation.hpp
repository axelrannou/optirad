#pragma once

#include <vector>

namespace optirad {

namespace Interpolation {
    double linear(double x0, double x1, double t);
    double bilinear(double v00, double v01, double v10, double v11, double tx, double ty);
    double trilinear(const double v[8], double tx, double ty, double tz);
}

} // namespace optirad
