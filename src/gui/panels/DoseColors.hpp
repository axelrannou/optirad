#pragma once

namespace optirad {

/// Shared color palette for dose statistics and DVH visualization.
struct DoseColors {
    static constexpr int kNumColors = 8;
    static constexpr float kColors[kNumColors][3] = {
        {1.0f, 0.0f, 0.0f},  // red
        {0.0f, 0.7f, 0.0f},  // green
        {0.2f, 0.4f, 1.0f},  // blue
        {1.0f, 0.8f, 0.0f},  // yellow
        {1.0f, 0.4f, 0.0f},  // orange
        {0.8f, 0.0f, 0.8f},  // purple
        {0.0f, 0.8f, 0.8f},  // cyan
        {0.6f, 0.3f, 0.1f},  // brown
    };
};

} // namespace optirad
