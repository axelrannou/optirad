#pragma once

#include "Patient.hpp"
#include "Volume.hpp"
#include "StructureSet.hpp"
#include <string>

namespace optirad {

class IDataImporter {
public:
    virtual ~IDataImporter() = default;

    virtual bool load(const std::string& path) = 0;

    virtual Patient getPatient() const = 0;
    virtual CTVolume getCT() const = 0;
    virtual StructureSet getStructureSet() const = 0;
};

} // namespace optirad
