# Utilities Module (`optirad_utils`)

The utilities module provides cross-cutting utility classes: logging, configuration, timing, math helpers, and interpolation.

**Library:** `optirad_utils`  
**Dependencies:** None (leaf library)

## Files

| File | Description |
|------|-------------|
| `Logger.hpp/cpp` | Thread-safe logging system |
| `Config.hpp/cpp` | JSON-based configuration |
| `Timer.hpp/cpp` | High-resolution timer |
| `MathUtils.hpp/cpp` | Vec3 math operations |
| `Interpolation.hpp/cpp` | Linear, bilinear, trilinear interpolation |

## Classes

### Logger

Thread-safe static logging system with severity levels.

```cpp
class Logger {
    static void init();
    static void info(const std::string& msg);
    static void warn(const std::string& msg);
    static void error(const std::string& msg);
    static void debug(const std::string& msg);     // Only in DEBUG builds
};
```

**Thread safety:** All methods are protected by a static `std::mutex`.

**Severity levels:**

| Level | Method | Description |
|-------|--------|-------------|
| INFO | `Logger::info()` | Normal operation messages |
| WARN | `Logger::warn()` | Non-critical issues |
| ERROR | `Logger::error()` | Errors that may affect results |
| DEBUG | `Logger::debug()` | Detailed diagnostics (debug builds only) |

**Usage:**
```cpp
Logger::init();
Logger::info("Loading DICOM from: " + path);
Logger::warn("No RT-STRUCT found, structures will be empty");
Logger::error("Failed to open file: " + filename);
```

### Config

JSON-based configuration loader/saver.

```cpp
class Config {
    bool load(const std::string& filepath);
    bool save(const std::string& filepath) const;

    std::string getDoseEngine() const;     // default: "PencilBeam"
    std::string getOptimizer() const;      // default: "LBFGS"
};
```

Used to persist and load default engine and optimizer selections.

### Timer

High-resolution timer for performance measurement.

```cpp
class Timer {
    Timer();                               // Auto-starts
    void start();
    void stop();
    double elapsedMs() const;              // Milliseconds
    double elapsedSeconds() const;
};
```

**Usage:**
```cpp
Timer timer;
timer.start();
// ... expensive computation ...
timer.stop();
Logger::info("Computation took " + std::to_string(timer.elapsedMs()) + " ms");
```

Uses `std::chrono::high_resolution_clock` internally.

### MathUtils (utils)

Basic 3D vector operations.

```cpp
using Vec3 = std::array<double, 3>;

namespace MathUtils {
    double dot(const Vec3& a, const Vec3& b);
    Vec3   cross(const Vec3& a, const Vec3& b);
    double norm(const Vec3& v);
    Vec3   normalize(const Vec3& v);
}
```

> **Note:** The geometry module provides a more comprehensive `MathUtils` in `geometry/MathUtils.hpp` with `Mat3`, rotation matrices, vector arithmetic, and matrix operations. The utils version provides the minimal subset for modules that don't depend on geometry.

### Interpolation

Standard interpolation functions.

```cpp
namespace Interpolation {
    // 1D linear interpolation: result = x0 + t × (x1 - x0)
    double linear(double x0, double x1, double t);

    // 2D bilinear interpolation over a unit square
    double bilinear(double v00, double v01, double v10, double v11,
                    double tx, double ty);

    // 3D trilinear interpolation over a unit cube
    double trilinear(const double v[8], double tx, double ty, double tz);
}
```

**Parameter conventions:**
- `t`, `tx`, `ty`, `tz` ∈ [0, 1] — normalized position within the cell
- `v[8]` — 8 corner values of the cube in binary order: `v[0]=v000, v[1]=v100, v[2]=v010, v[3]=v110, v[4]=v001, v[5]=v101, v[6]=v011, v[7]=v111`

**Usage:**
```cpp
// Trilinear interpolation within a voxel
double v[8] = { /* 8 corner values */ };
double result = Interpolation::trilinear(v, 0.3, 0.5, 0.7);
```

## Related Documentation

- [Geometry Module](geometry.md) — Extended MathUtils with Mat3 and rotation matrices
- [Dose Module](dose.md) — GridInterpolant2D for 2D kernel interpolation
