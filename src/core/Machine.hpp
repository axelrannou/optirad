#pragma once

#include <string>
#include <vector>
#include <array>
#include <cstdint>

namespace optirad {

struct MachineMeta {
    std::string radiationMode = "photons";
    std::string dataType = "-";
    std::string createdOn;
    std::string createdBy;
    std::string description;
    std::string name;
    double SAD = 1000.0;  // Source-Axis Distance (mm)
    double SCD = 500.0;   // Source-Collimator Distance (mm)
};

struct KernelEntry {
    double depth;
    std::vector<double> radialDistances;
    std::vector<double> primaryValues;
    std::vector<double> firstScatterValues;
    std::vector<double> secondScatterValues;
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

class Machine {
public:
    Machine() = default;

    static Machine createGenericPhoton();

    const MachineMeta& getMeta() const { return m_meta; }
    const MachineData& getData() const { return m_data; }
    const MachineConstraints& getConstraints() const { return m_constraints; }

    void setMeta(const MachineMeta& meta) { m_meta = meta; }
    void setData(const MachineData& data) { m_data = data; }
    void setConstraints(const MachineConstraints& constraints) { m_constraints = constraints; }

    const std::string& getName() const { return m_meta.name; }
    const std::string& getRadiationMode() const { return m_meta.radiationMode; }
    double getSAD() const { return m_meta.SAD; }
    double getSCD() const { return m_meta.SCD; }

private:
    MachineMeta m_meta;
    MachineData m_data;
    MachineConstraints m_constraints;
};

} // namespace optirad
