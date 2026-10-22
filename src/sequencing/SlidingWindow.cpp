#include "SlidingWindow.hpp"

#include <algorithm>
#include <cmath>
#include <utility>

namespace optirad {

namespace {

struct Exposed { int first = -1, last = -1; };

Exposed exposedRange(const std::vector<double>& profile) {
    Exposed range;
    for (size_t i = 0; i < profile.size(); ++i) {
        if (profile[i] <= 0.0) continue;
        if (range.first < 0) range.first = static_cast<int>(i);
        range.last = static_cast<int>(i);
    }
    return range;
}

using Path = std::vector<std::pair<double, double>>;   // (time, position) knots

double positionAt(const Path& path, double time) {
    if (time <= path.front().first) return path.front().second;
    if (time >= path.back().first) return path.back().second;
    const auto upper = std::upper_bound(path.begin(), path.end(), time,
        [](double t, const std::pair<double, double>& knot) { return t < knot.first; });
    const auto lower = upper - 1;
    const double span = upper->first - lower->first;
    if (span <= 0.0) return upper->second;
    return lower->second + (upper->second - lower->second) * (time - lower->first) / span;
}

} // namespace

double SlidingWindow::totalRise(const std::vector<double>& profile) {
    const Exposed range = exposedRange(profile);
    if (range.first < 0) return 0.0;
    double rises = profile[static_cast<size_t>(range.first)];
    for (int i = range.first; i < range.last; ++i)
        rises += std::max(0.0, profile[static_cast<size_t>(i) + 1] - profile[static_cast<size_t>(i)]);
    return rises;
}

double SlidingWindow::requiredMuPerSlot(const SlidingWindowRow& row, double speed, int slots) {
    const double rises = totalRise(row.profile);
    if (rises <= 0.0) return 0.0;
    const double sweepSlots = std::max(0.0, row.end - row.start) / speed;
    const double free = slots - sweepSlots;
    return free > 1e-9 ? rises / free : -1.0;
}

void SlidingWindow::sampleTips(const SlidingWindowRow& row, double speed, double muPerSlot, int slots,
                               std::vector<double>& trailing, std::vector<double>& leading) {
    trailing.assign(static_cast<size_t>(slots), row.start);
    leading.assign(static_cast<size_t>(slots), row.start);
    const Exposed range = exposedRange(row.profile);

    Path lead{{0.0, row.start}}, trail{{0.0, row.start}};
    double leadTime = 0.0, trailTime = 0.0;
    auto moveTo = [&](Path& path, double& time, double position) {
        time += std::max(0.0, position - path.back().second) / speed;
        path.emplace_back(time, position);
    };
    auto wait = [&](Path& path, double& time, double duration) {
        if (duration <= 0.0) return;
        time += duration;
        path.emplace_back(time, path.back().second);
    };

    if (range.first < 0 || muPerSlot <= 0.0) {
        moveTo(lead, leadTime, row.end);
        moveTo(trail, trailTime, row.end);
    } else {
        const auto& I = row.profile;
        const double edge0 = row.firstCenter + (range.first - 0.5) * row.pixel;
        moveTo(lead, leadTime, edge0);
        // Outside the exposed pixels the leaves travel together, closed.
        moveTo(trail, trailTime, edge0);
        wait(trail, trailTime, I[static_cast<size_t>(range.first)] / muPerSlot);
        for (int i = range.first; i <= range.last; ++i) {
            const double edge = row.firstCenter + (i + 0.5) * row.pixel;
            moveTo(lead, leadTime, edge);
            moveTo(trail, trailTime, edge);
            if (i == range.last) {
                wait(lead, leadTime, I[static_cast<size_t>(i)] / muPerSlot);
                break;
            }
            const double change = I[static_cast<size_t>(i) + 1] - I[static_cast<size_t>(i)];
            wait(lead, leadTime, std::max(0.0, -change) / muPerSlot);
            wait(trail, trailTime, std::max(0.0, change) / muPerSlot);
        }
        moveTo(lead, leadTime, row.end);
        moveTo(trail, trailTime, row.end);
    }

    for (int k = 0; k < slots; ++k) {
        const double centre = k + 0.5;
        trailing[static_cast<size_t>(k)] = positionAt(trail, centre);
        leading[static_cast<size_t>(k)] = positionAt(lead, centre);
    }
}

} // namespace optirad
