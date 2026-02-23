#pragma once

#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace optirad {

/// Machine type: Generic (pencil-beam kernel) vs PhaseSpace (IAEA PSF data)
enum class MachineType {
    Generic,      ///< Traditional pencil-beam kernel machine (e.g., machine_photons_Generic.json)
    PhaseSpace    ///< Phase-space based machine with IAEA PSF files (e.g., Varian TrueBeam)
};

struct MachineMeta {
    std::string radiationMode = "photons";
    std::string dataType = "-";
    std::string createdOn;
    std::string createdBy;
    std::string description;
    std::string name;
    double SAD = 1000.0;  // Source-Axis Distance (mm)
    double SCD = 500.0;   // Source-Collimator Distance (mm)
    MachineType machineType = MachineType::Generic;
};

struct KernelEntry {
    double SSD;                       // Source-Surface Distance (mm)
    std::vector<double> kernel1;      // Primary kernel values (one per kernelPos)
    std::vector<double> kernel2;      // First scatter kernel values
    std::vector<double> kernel3;      // Second scatter kernel values
};

struct MachineData {
    std::array<double, 3> betas = {0.3252, 0.0160, 0.0051};
    double energy = 6.0;  // MV
    double m = 0.0051;
    std::vector<std::array<double, 2>> primaryFluence;
    std::vector<KernelEntry> kernel;
    std::vector<double> kernelPos;
    double penumbraFWHMatIso = 5.0;  // mm
};

struct MachineConstraints {
    std::array<double, 2> gantryRotationSpeed = {0.0, 6.0};   // deg/s
    std::array<double, 2> leafSpeed = {0.0, 60.0};             // mm/s
    std::array<double, 2> monitorUnitRate = {1.25, 10.0};      // MU/s
};

/// Geometry description for machines with jaws, MLC, and collimator.
/// Populated from folder-based machine JSON (e.g., Varian_TrueBeam6MV.json).
struct MachineGeometry {
    // Jaw limits (mm)
    double jawX1Min = 0.0, jawX1Max = 200.0;
    double jawX2Min = 0.0, jawX2Max = 200.0;
    double jawY1Min = 0.0, jawY1Max = 200.0;
    double jawY2Min = 0.0, jawY2Max = 200.0;
    std::array<double, 2> defaultFieldSize = {100.0, 100.0}; // mm

    // Collimator angle range (deg)
    std::array<double, 2> collimatorRange = {-180.0, 180.0};
    double defaultCollimatorAngle = 0.0;

    // Couch angle range (deg)
    std::array<double, 2> couchRange = {-90.0, 90.0};
    double defaultCouchAngle = 0.0;

    // MLC
    std::string mlcType;
    int numLeaves = 0;
    std::vector<double> leafWidths; // mm (e.g., 5.0 and 10.0)
    double maxLeafTravel = 0.0;     // mm
    bool interdigitation = false;

    // Phase-space files directory (absolute path to folder containing .IAEAphsp/.IAEAheader pairs)
    std::string phaseSpaceDir;
    int numPhaseSpaceFiles = 0;
    std::vector<std::string> phaseSpaceFileNames; // base names without extension

    // Beam energy
    double beamEnergyMV = 6.0;
    double doseRateMUPerMin = 600.0;
    double apertureSamplingResolutionMm = 2.0;
};

class Machine {
public:
    Machine() = default;

    const MachineMeta& getMeta() const { return m_meta; }
    const MachineData& getData() const { return m_data; }
    const MachineConstraints& getConstraints() const { return m_constraints; }
    const MachineGeometry& getGeometry() const { return m_geometry; }

    void setMeta(const MachineMeta& meta) { m_meta = meta; }
    void setData(const MachineData& data) { m_data = data; }
    void setConstraints(const MachineConstraints& constraints) { m_constraints = constraints; }
    void setGeometry(const MachineGeometry& geometry) { m_geometry = geometry; }

    const std::string& getName() const { return m_meta.name; }
    const std::string& getRadiationMode() const { return m_meta.radiationMode; }
    double getSAD() const { return m_meta.SAD; }
    double getSCD() const { return m_meta.SCD; }
    MachineType getMachineType() const { return m_meta.machineType; }
    bool isPhaseSpace() const { return m_meta.machineType == MachineType::PhaseSpace; }

private:
    MachineMeta m_meta;
    MachineData m_data;
    MachineConstraints m_constraints;
    MachineGeometry m_geometry;
};

} // namespace optirad
