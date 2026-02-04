#pragma once

#include <string>
#include <memory>

namespace optirad {

class Patient;
class StructureSet;

class IDataImporter {
public:
    virtual ~IDataImporter() = default;

    virtual bool canImport(const std::string& path) const = 0;
    virtual std::unique_ptr<Patient> importPatient(const std::string& path) = 0;
    virtual std::unique_ptr<StructureSet> importStructures(const std::string& path) = 0;
};

} // namespace optirad
