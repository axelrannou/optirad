#include "Interpolation.hpp"

namespace optirad::Interpolation {

double linear(double x0, double x1, double t) {
    return x0 * (1.0 - t) + x1 * t;
}

double bilinear(double v00, double v01, double v10, double v11, double tx, double ty) {
    return linear(linear(v00, v01, tx), linear(v10, v11, tx), ty);
}

double trilinear(const double v[8], double tx, double ty, double tz) {
    double c0 = bilinear(v[0], v[1], v[2], v[3], tx, ty);
    double c1 = bilinear(v[4], v[5], v[6], v[7], tx, ty);
    return linear(c0, c1, tz);
}

} // namespace optirad::Interpolation
