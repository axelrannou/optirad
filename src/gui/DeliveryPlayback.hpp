#pragma once

#include "core/Aperture.hpp"
#include "core/ArcPlan.hpp"
#include "core/Beam.hpp"
#include <memory>
#include <string>
#include <vector>

namespace optirad {

struct GuiAppState;
struct ArcPipelineResult;
class Plan;
class Stf;
class PatientData;

struct DeliveryLine {
    Vec3 from;
    Vec3 to;
};

struct DeliveryOverlay {
    std::vector<DeliveryLine> path;
    std::vector<DeliveryLine> leaves;
    std::vector<DeliveryLine> aim;
    std::vector<Vec3> fitPoints;
};

struct DeliveryFrame {
    bool valid = false;
    bool beamOn = true;
    bool hasAperture = false;
    size_t group = 0;
    size_t beamIndex = 0;
    size_t controlPoint = 0;
    int segment = -1;
    Beam beam;
    Aperture aperture;
    std::vector<double> boundaries;
    double gantrySpeed = 0.0;
    double doseRate = 0.0; // machine MU/fraction/s
    double leafSpeed = 0.0;
};

/// One plan-wide clock. Only Application advances it; views only seek/read it.
class DeliveryPlayback {
public:
    void synchronize(const GuiAppState& state);
    void advance(double seconds);
    void seek(double seconds);
    void rewind();
    void togglePlaying();
    void selectGroup(size_t group);
    void selectSegment(size_t beam, int segment);

    bool available() const { return !m_steps.empty(); }
    bool isVmat() const { return m_vmat; }
    bool estimated() const { return m_estimated; }
    double time() const { return m_time; }
    double duration() const { return m_duration; }
    size_t groupCount() const { return m_groupStarts.size(); }
    const DeliveryFrame& frame() const { return m_frame; }
    const DeliveryOverlay& overlay() const { return m_overlay; }
    const std::string& description() const { return m_description; }
    double extentX() const { return m_extentX; }
    double extentZ() const { return m_extentZ; }

    bool playing = false;
    bool showOverlays = true;
    float playbackSpeed = 1.0f;

private:
    struct Step {
        DeliveryFrame frame;
        Aperture before;
        Aperture after;
        double start = 0.0;
        double duration = 0.0;
        double angleStart = 0.0;
        double angleDelta = 0.0;
        bool movingLeaves = false;
    };
    struct Path {
        ArcSpec arc;
        Beam beam;
        std::vector<double> controlPointAngles;
    };

    void rebuild(const GuiAppState& state);
    void sample();
    void buildPaths();
    void buildOverlay();
    void append(Step step);

    std::vector<Step> m_steps;
    std::vector<Path> m_paths;
    std::vector<double> m_groupStarts;
    std::vector<DeliveryLine> m_pathLines;
    std::vector<Vec3> m_fitPoints;
    DeliveryFrame m_frame;
    DeliveryOverlay m_overlay;
    double m_time = 0.0, m_duration = 0.0;
    double m_extentX = 20.0, m_extentZ = 20.0;
    bool m_vmat = false, m_estimated = true;
    std::string m_description;

    std::shared_ptr<Plan> m_plan;
    std::shared_ptr<Stf> m_stf;
    std::shared_ptr<PatientData> m_patient;
    std::shared_ptr<ArcPipelineResult> m_result;
    int m_doseId = -1;
    int m_doseVersion = -1;
    size_t m_resultTiming = 0;
    size_t m_resultCount = 0;
};

} // namespace optirad
