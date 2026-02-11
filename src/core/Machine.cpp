#include "Machine.hpp"
#include <cmath>

namespace optirad {

Machine Machine::createGenericPhoton() {
    Machine machine;

    // Meta
    MachineMeta meta;
    meta.radiationMode = "photons";
    meta.dataType = "-";
    meta.createdOn = "27-Oct-2015";
    meta.createdBy = "wieserh";
    meta.description = "photon pencil beam kernels for a 6MV machine";
    meta.name = "Generic";
    meta.SAD = 1000.0;
    meta.SCD = 500.0;
    machine.setMeta(meta);

    // Data
    MachineData data;
    data.betas = {0.3252, 0.0160, 0.0051};
    data.energy = 6.0;
    data.m = 0.0051;
    data.penumbraFWHMatIso = 5.0;

    // Generate kernel positions: 0.0, 0.5, 1.0, ..., 179.5 (360 entries)
    data.kernelPos.reserve(360);
    for (int i = 0; i < 360; ++i) {
        data.kernelPos.push_back(i * 0.5);
    }

    // Primary fluence placeholder (38 entries as in matRad)
    // Simplified: linear ramp for off-axis ratio
    data.primaryFluence.reserve(38);
    for (int i = 0; i < 38; ++i) {
        double offAxis = i * 10.0; // mm
        double fluence = std::exp(-0.0001 * offAxis * offAxis);
        data.primaryFluence.push_back({offAxis, fluence});
    }

    machine.setData(data);

    // Constraints
    MachineConstraints constraints;
    constraints.gantryRotationSpeed = {0.0, 6.0};
    constraints.leafSpeed = {0.0, 60.0};
    constraints.monitorUnitRate = {1.25, 10.0};
    machine.setConstraints(constraints);

    return machine;
}

} // namespace optirad
