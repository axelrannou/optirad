#pragma once

#include "Machine.hpp"
#include "../steering/StfProperties.hpp"
#include "Beam.hpp"
#include "PatientData.hpp"

#include <string>
#include <vector>
#include <memory>
#include <array>

namespace optirad {

class Plan {
public:
    Plan() = default;

    // Name
    void setName(const std::string& name);
    const std::string& getName() const;

    // Radiation mode
    void setRadiationMode(const std::string& mode) { m_radiationMode = mode; }
    const std::string& getRadiationMode() const { return m_radiationMode; }

    // Machine
    void setMachine(const Machine& machine) { m_machine = machine; }
    const Machine& getMachine() const { return m_machine; }

    // Fractions
    void setNumOfFractions(int n) { m_numOfFractions = n; }
    int getNumOfFractions() const { return m_numOfFractions; }

    // STF properties
    void setStfProperties(const StfProperties& stf) { m_stfProperties = stf; }
    const StfProperties& getStfProperties() const { return m_stfProperties; }
    StfProperties& stfProperties() { return m_stfProperties; }

    // Patient data reference
    void setPatientData(std::shared_ptr<PatientData> data) { m_patientData = data; }
    std::shared_ptr<PatientData> getPatientData() const { return m_patientData; }

    // Beams (kept for compatibility)
    void addBeam(const Beam& beam);
    const std::vector<Beam>& getBeams() const;
    size_t getNumBeams() const;

    /// Compute isocenter from CT volume center of mass (simplified)
    std::array<double, 3> computeIsoCenter() const;

    // Dose grid resolution
    void setDoseGridResolution(const std::array<double, 3>& res) { m_doseGridResolution = res; }
    const std::array<double, 3>& getDoseGridResolution() const { return m_doseGridResolution; }

    // Prescribed dose (single source of truth)
    void setPrescribedDose(double dose) { m_prescribedDose = dose; }
    double getPrescribedDose() const { return m_prescribedDose; }

    /// Print plan summary
    void printSummary() const;

private:
    std::string m_name;
    std::string m_radiationMode = "photons";
    Machine m_machine;
    int m_numOfFractions = 30;
    StfProperties m_stfProperties;
    std::shared_ptr<PatientData> m_patientData;
    std::vector<Beam> m_beams;
    std::array<double, 3> m_doseGridResolution = {2.5, 2.5, 2.5}; // mm
    double m_prescribedDose = 0.0; // Gy (0 = not set, use default)
};

} // namespace optirad
