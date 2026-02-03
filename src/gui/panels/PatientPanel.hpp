#pragma once

#include "IPanel.hpp"
#include "core/Patient.hpp"
#include "geometry/StructureSet.hpp"

namespace optirad {

class PatientPanel : public IPanel {
public:
    std::string getName() const override { return "Patient"; }
    void update() override;
    void render() override;

    void setPatient(const Patient* patient) { m_patient = patient; }
    void setStructureSet(const StructureSet* structureSet) { m_structureSet = structureSet; }

private:
    const Patient* m_patient = nullptr;
    const StructureSet* m_structureSet = nullptr;
};

} // namespace optirad
