#include "DeliveryPlayback.hpp"
#include "AppState.hpp"
#include "utils/Logger.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <stdexcept>

namespace optirad {
namespace {

double wrap(double angle) {
    angle = std::fmod(angle, 360.0);
    return angle < 0.0 ? angle + 360.0 : angle;
}

Beam geometry(const Beam& source) {
    Beam beam;
    beam.setIsocenter(source.getIsocenter());
    beam.setGantryAngle(source.getGantryAngle());
    beam.setCouchAngle(source.getCouchAngle());
    beam.setCollimatorAngle(source.getCollimatorAngle());
    beam.setSAD(source.getSAD());
    beam.setSCD(source.getSCD());
    beam.setFieldSize(source.getFieldSize());
    beam.setBixelWidth(source.getBixelWidth());
    return beam;
}

Vec3 point(const Beam& beam, double x, double y, double z) {
    return vecAdd(beam.getIsocenter(),
                  beam.getBeamRotationMatrix() * Vec3{x, y, z});
}

bool matches(const Aperture& a, const Aperture& b) {
    return a.bankA.size() == a.bankB.size() && b.bankA.size() == b.bankB.size() &&
           a.bankA.size() == b.bankA.size();
}

void interpolate(Aperture& result, const Aperture& from, const Aperture& to, double t) {
    for (size_t pair = 0; pair < result.bankA.size(); ++pair) {
        result.bankA[pair] = from.bankA[pair] + t * (to.bankA[pair] - from.bankA[pair]);
        result.bankB[pair] = from.bankB[pair] + t * (to.bankB[pair] - from.bankB[pair]);
    }
}

} // namespace

void DeliveryPlayback::synchronize(const GuiAppState& state) {
    const auto* dose = state.doseStore.getSelected();
    const int id = dose ? dose->id : -1;
    const auto cached = state.arcCache.find(id);
    const auto result = cached != state.arcCache.end() ? cached->second.result : nullptr;
    size_t resultTiming = result ? std::hash<double>{}(result->sequence.totalTimeSec) : 0;
    if (result) {
        for (const auto& cp : result->sequence.controlPoints)
            resultTiming ^= std::hash<double>{}(cp.timeSec) + 0x9e3779b9 +
                            (resultTiming << 6) + (resultTiming >> 2);
    }
    const size_t resultCount = result ? result->sequence.controlPoints.size() : 0;
    if (m_plan == state.plan && m_stf == state.stf && m_patient == state.patientData &&
        m_result == result && m_doseId == id && m_doseVersion == state.doseStore.version() &&
        m_resultTiming == resultTiming && m_resultCount == resultCount) return;
    m_plan = state.plan;
    m_stf = state.stf;
    m_patient = state.patientData;
    m_result = result;
    m_doseId = id;
    m_doseVersion = state.doseStore.version();
    m_resultTiming = resultTiming;
    m_resultCount = resultCount;
    try {
        rebuild(state);
    } catch (const std::invalid_argument& error) {
        m_description = std::string("Cannot preview delivery: ") + error.what();
        Logger::error(m_description);
        m_steps.clear();
        m_groupStarts.clear();
        m_frame = {};
        m_overlay = {};
        m_duration = m_time = 0.0;
        playing = false;
    }
}

void DeliveryPlayback::append(Step step) {
    if (!std::isfinite(step.duration) || step.duration < 0.0)
        throw std::invalid_argument("Delivery timeline has invalid control-point timing");
    const auto& beam = step.frame.beam;
    if (!std::isfinite(beam.getSAD()) || !std::isfinite(beam.getSCD()) ||
        beam.getSCD() <= 0.0 || beam.getSAD() <= beam.getSCD())
        throw std::invalid_argument("beam geometry requires 0 < SCD < SAD");
    if (!std::isfinite(step.angleStart) || !std::isfinite(step.angleDelta) ||
        !std::isfinite(beam.getCouchAngle()) ||
        std::any_of(beam.getIsocenter().begin(), beam.getIsocenter().end(),
                    [](double coordinate) { return !std::isfinite(coordinate); }))
        throw std::invalid_argument("Delivery timeline has invalid beam coordinates");
    step.start = m_duration;
    m_duration += step.duration;
    m_extentX = std::max(m_extentX, step.frame.beam.getFieldSize()[0] * 0.5 + 15.0);
    m_extentZ = std::max(m_extentZ, step.frame.beam.getFieldSize()[1] * 0.5 + 5.0);
    if (step.frame.hasAperture || !step.frame.aperture.bankB.empty()) {
        const auto& aperture = step.frame.aperture;
        if (aperture.bankA.size() != aperture.bankB.size() ||
            step.frame.boundaries.size() != aperture.bankA.size() + 1)
            throw std::invalid_argument("Delivery timeline has inconsistent MLC geometry");
        for (size_t pair = 0; pair < aperture.bankA.size(); ++pair) {
            if (!std::isfinite(aperture.bankA[pair]) || !std::isfinite(aperture.bankB[pair]) ||
                !std::isfinite(step.frame.boundaries[pair]) ||
                !std::isfinite(step.frame.boundaries[pair + 1]) ||
                step.frame.boundaries[pair + 1] <= step.frame.boundaries[pair])
                throw std::invalid_argument("Delivery timeline has invalid MLC coordinates");
        }
        for (double tip : aperture.bankA) m_extentX = std::max(m_extentX, std::abs(tip) + 15.0);
        for (double tip : aperture.bankB) m_extentX = std::max(m_extentX, std::abs(tip) + 15.0);
        for (double z : step.frame.boundaries) m_extentZ = std::max(m_extentZ, std::abs(z) + 5.0);
    }
    if (step.duration > 0.0) m_steps.push_back(std::move(step));
}

void DeliveryPlayback::rebuild(const GuiAppState& state) {
    playing = false;
    m_time = m_duration = 0.0;
    m_extentX = m_extentZ = 20.0;
    m_steps.clear();
    m_paths.clear();
    m_groupStarts.clear();
    m_pathLines.clear();
    m_fitPoints.clear();
    m_frame = {};
    m_overlay = {};
    m_vmat = m_result || (m_plan && m_plan->isVmat());
    m_estimated = !m_result || m_result->sequence.controlPoints.empty();
    m_description.clear();
    const MachineConstraints limits = m_plan ? m_plan->getMachine().getConstraints() : MachineConstraints{};
    const double gantryMax = limits.gantryRotationSpeed[1];
    if (!std::isfinite(gantryMax) || gantryMax <= 0.0) {
        m_description = "Cannot preview delivery: invalid machine gantry speed.";
        Logger::error(m_description);
        return;
    }
    if (m_result && !m_result->sequence.controlPoints.empty()) {
        const auto& points = m_result->sequence.controlPoints;
        if (!m_result->controlPointStf || m_result->controlPointStf->getCount() != points.size()) {
            m_description = "Cannot preview delivery: missing control-point beam geometry.";
            Logger::error(m_description);
            return;
        }
        m_description = "Planned VMAT timing; inter-arc repositioning is not modeled.";
        const double muScale = m_plan ? m_plan->machineMuPerWeight() : 1.0;
        for (size_t begin = 0; begin < points.size();) {
            size_t end = begin + 1;
            while (end < points.size() && points[end].geometry.arcIndex == points[begin].geometry.arcIndex) ++end;
            bool clockwise = true;
            if (m_plan && points[begin].geometry.arcIndex < m_plan->getVmatSettings().arcs.size())
                clockwise = m_plan->getVmatSettings().arcs[points[begin].geometry.arcIndex].clockwise;
            if (end > begin + 1) {
                const double delta = std::remainder(points[begin + 1].geometry.gantryAngle -
                                                     points[begin].geometry.gantryAngle, 360.0);
                if (std::abs(delta) > 1e-8 && std::abs(delta) < 180.0) clockwise = delta > 0.0;
            }
            const double sign = clockwise ? 1.0 : -1.0;
            Path path;
            path.beam = geometry(*m_result->controlPointStf->getBeam(begin));
            path.arc.startAngle = wrap(points[begin].geometry.gantryAngle - sign * points[begin].geometry.width * 0.5);
            path.arc.couchAngle = points[begin].geometry.couchAngle;
            path.arc.clockwise = clockwise;
            double span = 0.0;
            m_groupStarts.push_back(m_duration);
            for (size_t i = begin; i < end; ++i) {
                const auto& cp = points[i];
                Step step;
                step.frame.valid = true;
                step.frame.group = m_paths.size();
                step.frame.beamIndex = cp.geometry.sectorIndex;
                step.frame.controlPoint = i;
                step.frame.beam = geometry(*m_result->controlPointStf->getBeam(i));
                step.frame.aperture = cp.aperture;
                step.frame.beamOn = cp.mu > 0.0;
                step.frame.hasAperture = !cp.aperture.bankA.empty();
                step.frame.boundaries = cp.leafPairBoundariesZ;
                step.frame.gantrySpeed = cp.timeSec > 0.0 ? cp.geometry.width / cp.timeSec : 0.0;
                step.frame.doseRate = cp.doseRate * muScale;
                step.frame.leafSpeed = cp.leafSpeed;
                step.duration = cp.timeSec;
                step.angleStart = cp.geometry.gantryAngle - sign * cp.geometry.width * 0.5;
                step.angleDelta = sign * cp.geometry.width;
                step.before = step.after = cp.aperture;
                if (i > begin && points[i - 1].leafPairBoundariesZ == cp.leafPairBoundariesZ &&
                    matches(points[i - 1].aperture, cp.aperture))
                    interpolate(step.before, points[i - 1].aperture, cp.aperture, 0.5);
                if (i + 1 < end && points[i + 1].leafPairBoundariesZ == cp.leafPairBoundariesZ &&
                    matches(points[i + 1].aperture, cp.aperture))
                    interpolate(step.after, cp.aperture, points[i + 1].aperture, 0.5);
                step.movingLeaves = step.frame.hasAperture;
                append(std::move(step));
                path.controlPointAngles.push_back(cp.geometry.gantryAngle);
                span += cp.geometry.width;
            }
            path.arc.stopAngle = wrap(path.arc.startAngle + sign * span);
            m_paths.push_back(std::move(path));
            begin = end;
        }
    } else if (m_vmat && m_plan && m_stf && !m_stf->isEmpty()) {
        m_description = "Geometry preview at the gantry speed limit (no delivered MLC results).";
        const auto& arcs = m_plan->getVmatSettings().arcs;
        const auto sectors = ArcGeometry::sectors(m_plan->getVmatSettings());
        for (size_t arc = 0; arc < arcs.size(); ++arc) {
            auto first = std::find_if(sectors.begin(), sectors.end(),
                [arc](const ArcSector& sector) { return sector.arcIndex == arc; });
            if (first == sectors.end()) continue;
            const auto beamIndex = static_cast<size_t>(first - sectors.begin());
            const auto* source = m_stf->getBeam(beamIndex);
            if (!source) throw std::invalid_argument("missing beam geometry for a planned arc");
            Step step;
            step.frame.valid = true;
            step.frame.group = m_paths.size();
            step.frame.beamIndex = beamIndex;
            step.frame.beam = geometry(*source);
            step.frame.gantrySpeed = gantryMax;
            step.duration = arcs[arc].span() / gantryMax;
            step.angleStart = arcs[arc].startAngle;
            step.angleDelta = (arcs[arc].clockwise ? 1.0 : -1.0) * arcs[arc].span();
            m_groupStarts.push_back(m_duration);
            Path path{arcs[arc], step.frame.beam, {}};
            for (const auto& sector : sectors)
                if (sector.arcIndex == arc) path.controlPointAngles.push_back(sector.centerAngle);
            m_paths.push_back(std::move(path));
            append(std::move(step));
        }
    } else if (m_stf && !m_stf->isEmpty()) {
        const auto cache = state.seqCache.find(m_doseId);
        const auto* sequences = cache != state.seqCache.end() ? &cache->second.sequences
                              : m_plan && m_plan == state.importedPlan && !state.importedSequences.empty()
                                  ? &state.importedSequences : nullptr;
        m_description = sequences ? "Estimated SMLC timing; transitions are beam-off, couch motion is not modeled."
                                  : "Geometry preview without MLC results: illustrative 1 s dwell per beam.";
        const double muScale = m_plan ? m_plan->machineMuPerWeight() : 1.0;
        const double rate = limits.monitorUnitRate[1] / muScale;
        if (!std::isfinite(rate) || rate <= 0.0 ||
            !std::isfinite(limits.leafSpeed[1]) || limits.leafSpeed[1] <= 0.0) {
            m_description = "Cannot preview delivery: invalid machine dose rate or leaf speed.";
            Logger::error(m_description);
            return;
        }
        for (size_t beamIndex = 0; beamIndex < m_stf->getCount(); ++beamIndex) {
            const Beam beam = geometry(*m_stf->getBeam(beamIndex));
            m_groupStarts.push_back(m_duration);
            const auto* sequence = sequences && beamIndex < sequences->size() ? &(*sequences)[beamIndex] : nullptr;
            const size_t count = sequence && !sequence->segments.empty() ? sequence->segments.size() : 1;
            for (size_t segment = 0; segment < count; ++segment) {
                Step step;
                step.frame.valid = true;
                step.frame.group = beamIndex;
                step.frame.beamIndex = beamIndex;
                step.frame.beam = beam;
                step.frame.hasAperture = sequence && !sequence->segments.empty();
                step.frame.segment = step.frame.hasAperture ? static_cast<int>(segment) : -1;
                step.angleStart = beam.getGantryAngle();
                if (step.frame.hasAperture) {
                    step.frame.aperture = sequence->segments[segment];
                    step.frame.boundaries = sequence->leafPairBoundariesZ;
                    if (step.frame.boundaries.empty()) {
                        double firstZ = 0.0;
                        if (!m_stf->getBeam(beamIndex)->getRays().empty()) {
                            firstZ = m_stf->getBeam(beamIndex)->getRays().front().getRayPosBev()[2];
                            for (const auto& ray : m_stf->getBeam(beamIndex)->getRays())
                                firstZ = std::min(firstZ, ray.getRayPosBev()[2]);
                        }
                        for (size_t pair = 0; pair <= step.frame.aperture.bankA.size(); ++pair)
                            step.frame.boundaries.push_back(firstZ + (static_cast<double>(pair) - 0.5) * beam.getBixelWidth());
                    }
                    step.duration = step.frame.aperture.weight / rate;
                    step.frame.doseRate = limits.monitorUnitRate[1];
                } else {
                    step.duration = 1.0;
                }
                if (!m_steps.empty()) {
                    const auto& previous = m_steps.back().frame;
                    Step transit = step;
                    transit.frame.beamOn = false;
                    transit.frame.doseRate = 0.0;
                    transit.angleStart = previous.beam.getGantryAngle();
                    transit.angleDelta = std::remainder(beam.getGantryAngle() - transit.angleStart, 360.0);
                    transit.duration = std::abs(transit.angleDelta) / gantryMax;
                    if (step.frame.hasAperture && previous.boundaries == step.frame.boundaries &&
                        matches(previous.aperture, step.frame.aperture)) {
                        double travel = 0.0;
                        for (size_t pair = 0; pair < step.frame.aperture.bankA.size(); ++pair)
                            travel = std::max({travel,
                                std::abs(previous.aperture.bankA[pair] - step.frame.aperture.bankA[pair]),
                                std::abs(previous.aperture.bankB[pair] - step.frame.aperture.bankB[pair])});
                        transit.duration = std::max(transit.duration, travel / limits.leafSpeed[1]);
                        transit.before = previous.aperture;
                        transit.after = step.frame.aperture;
                        transit.movingLeaves = true;
                        transit.frame.leafSpeed = transit.duration > 0.0 ? travel / transit.duration : 0.0;
                    }
                    transit.frame.gantrySpeed = transit.duration > 0.0 ?
                        std::abs(transit.angleDelta) / transit.duration : 0.0;
                    append(std::move(transit));
                }
                append(std::move(step));
            }
        }
    }
    if (!available()) {
        if (m_description.empty()) m_description = "Create a plan to preview delivery.";
        else m_description += " No positive-time delivery steps are available.";
    }
    buildPaths();
    sample();
}

void DeliveryPlayback::seek(double seconds) {
    if (!std::isfinite(seconds)) throw std::invalid_argument("Delivery seek time must be finite");
    m_time = std::clamp(seconds, 0.0, m_duration);
    if (m_time >= m_duration) playing = false;
    sample();
}

void DeliveryPlayback::advance(double seconds) {
    if (!std::isfinite(seconds) || seconds < 0.0)
        throw std::invalid_argument("Delivery frame duration must be finite and nonnegative");
    if (!std::isfinite(playbackSpeed) || playbackSpeed <= 0.0f)
        throw std::invalid_argument("Delivery playback multiplier must be finite and positive");
    if (playing && available()) seek(m_time + seconds * playbackSpeed);
}

void DeliveryPlayback::rewind() {
    playing = false;
    seek(0.0);
}

void DeliveryPlayback::togglePlaying() {
    if (!available()) return;
    if (!playing && m_time >= m_duration) seek(0.0);
    playing = !playing;
}

void DeliveryPlayback::selectGroup(size_t group) {
    if (group >= m_groupStarts.size()) throw std::out_of_range("Delivery group is out of range");
    playing = false;
    seek(m_groupStarts[group]);
}

void DeliveryPlayback::selectSegment(size_t beam, int segment) {
    const auto found = std::find_if(m_steps.begin(), m_steps.end(), [beam, segment](const Step& step) {
        return step.frame.beamIndex == beam && step.frame.segment == segment && step.frame.beamOn;
    });
    if (found == m_steps.end()) throw std::out_of_range("Delivery segment is out of range");
    playing = false;
    seek(found->start);
}

void DeliveryPlayback::sample() {
    if (m_steps.empty()) return;
    auto step = std::upper_bound(m_steps.begin(), m_steps.end(), m_time,
        [](double time, const Step& item) { return time < item.start; });
    if (step != m_steps.begin()) --step;
    const double t = std::clamp((m_time - step->start) / step->duration, 0.0, 1.0);
    m_frame = step->frame;
    m_frame.beam.setGantryAngle(wrap(step->angleStart + t * step->angleDelta));
    if (step->movingLeaves) {
        if (!m_vmat) interpolate(m_frame.aperture, step->before, step->after, t);
        else if (t <= 0.5) interpolate(m_frame.aperture, step->before, step->frame.aperture, t * 2.0);
        else interpolate(m_frame.aperture, step->frame.aperture, step->after, (t - 0.5) * 2.0);
    }
    buildOverlay();
}

void DeliveryPlayback::buildPaths() {
    for (const auto& path : m_paths) {
        Beam beam = path.beam;
        beam.setCouchAngle(path.arc.couchAngle);
        const double radius = beam.getSAD() - beam.getSCD();
        const int steps = std::max(2, static_cast<int>(std::ceil(path.arc.span() / 2.0)));
        Vec3 previous{};
        for (int i = 0; i <= steps; ++i) {
            beam.setGantryAngle(path.arc.angleAt(path.arc.span() * i / steps));
            const auto center = point(beam, 0.0, -radius, 0.0);
            if (i > 0) m_pathLines.push_back({previous, center});
            previous = center;
            const double scale = beam.getSCD() / beam.getSAD();
            for (double x : {-m_extentX, m_extentX})
                for (double z : {-m_extentZ, m_extentZ})
                    m_fitPoints.push_back(point(beam, x * scale, -radius, z * scale));
        }
        const auto tip = previous;
        const auto tangent = normalize(vecSub(tip, m_pathLines.back().from));
        const auto inward = normalize(vecSub(beam.getIsocenter(), tip));
        m_pathLines.push_back({tip, vecAdd(vecSub(tip, vecScale(tangent, 25.0)), vecScale(inward, 12.0))});
        m_pathLines.push_back({tip, vecSub(vecSub(tip, vecScale(tangent, 25.0)), vecScale(inward, 12.0))});
        // Radial tick per control point, straddling the arc in the gantry plane.
        for (double angle : path.controlPointAngles) {
            beam.setGantryAngle(angle);
            m_pathLines.push_back({point(beam, 0.0, -radius - 8.0, 0.0), point(beam, 0.0, -radius + 8.0, 0.0)});
        }
    }
    if (m_paths.empty()) {
        for (const auto& step : m_steps) {
            auto beam = step.frame.beam;
            const double scale = beam.getSCD() / beam.getSAD();
            const int count = std::max(1, static_cast<int>(std::ceil(std::abs(step.angleDelta) / 2.0)));
            for (int i = 0; i <= count; ++i) {
                beam.setGantryAngle(step.angleStart + step.angleDelta * i / count);
                for (double x : {-m_extentX, m_extentX})
                    for (double z : {-m_extentZ, m_extentZ})
                        m_fitPoints.push_back(point(beam, x * scale, beam.getSCD() - beam.getSAD(), z * scale));
            }
        }
    }
}

void DeliveryPlayback::buildOverlay() {
    m_overlay.path = m_pathLines;
    m_overlay.fitPoints = m_fitPoints;
    m_overlay.leaves.clear();
    m_overlay.aim.clear();
    const auto& beam = m_frame.beam;
    const double scale = beam.getSCD() / beam.getSAD();
    auto leaf = [&](double x, double z) {
        return point(beam, x * scale, beam.getSCD() - beam.getSAD(), z * scale);
    };
    auto rectangle = [&](double x0, double x1, double z0, double z1) {
        const auto a = leaf(x0, z0), b = leaf(x1, z0), c = leaf(x1, z1), d = leaf(x0, z1);
        m_overlay.leaves.insert(m_overlay.leaves.end(), {{a, b}, {b, c}, {c, d}, {d, a}});
    };
    auto aim = [&](double x, double z) {
        const double depth = 150.0;
        const double divergence = (beam.getSAD() + depth) / beam.getSAD();
        m_overlay.aim.push_back({leaf(x, z), point(beam, x * divergence, depth, z * divergence)});
    };
    if (m_frame.hasAperture) {
        const auto& aperture = m_frame.aperture;
        const size_t count = aperture.bankA.size();
        for (size_t pair = 0; pair < count; ++pair) {
            const double a = aperture.bankA[pair], b = aperture.bankB[pair];
            const double z0 = m_frame.boundaries[pair], z1 = m_frame.boundaries[pair + 1];
            rectangle(-m_extentX, a, z0, z1);
            rectangle(b, m_extentX, z0, z1);
            if (b > a && pair % std::max(size_t{1}, count / 4) == 0) {
                aim(a, (z0 + z1) * 0.5);
                aim(b, (z0 + z1) * 0.5);
            }
        }
    } else {
        const auto field = beam.getFieldSize();
        rectangle(-field[0] * 0.5, field[0] * 0.5, -field[1] * 0.5, field[1] * 0.5);
        for (double x : {-field[0] * 0.5, field[0] * 0.5})
            for (double z : {-field[1] * 0.5, field[1] * 0.5}) aim(x, z);
    }
    m_overlay.aim.push_back({leaf(0.0, 0.0), point(beam, 0.0, 150.0, 0.0)});
    const auto iso = beam.getIsocenter();
    for (size_t axis = 0; axis < 3; ++axis) {
        Vec3 a = iso, b = iso;
        a[axis] -= 6.0; b[axis] += 6.0;
        m_overlay.aim.push_back({a, b});
    }
}

} // namespace optirad
