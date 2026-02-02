#pragma once

#include "IPanel.hpp"
#include "Patient.hpp"
#include "StructureSet.hpp"

namespace optirad {

class PatientPanel : public IPanel {
public:
    std::string getName() const override { return "Patient"; }
    void update() override;
    void render() override;

    void setPatient(const Patient& patient);
    void setStructureSet(const StructureSet& structureSet);

private:
    Patient m_patient;
    StructureSet m_structureSet;
};

} // namespace optirad
