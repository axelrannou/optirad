#pragma once

#include "IDataImporter.hpp"

namespace optirad {

class DicomImporter : public IDataImporter {
public:
    bool load(const std::string& path) override;

    Patient getPatient() const override;
    CTVolume getCT() const override;
    StructureSet getStructureSet() const override;

private:
    Patient m_patient;
    CTVolume m_ct;
    StructureSet m_structureSet;
};

} // namespace optirad
