#pragma once

#include "core/Machine.hpp"
#include <nlohmann/json_fwd.hpp>
#include <string>
#include <vector>

namespace optirad {

/// Loads machine data from JSON files.
///
/// Two path conventions are supported:
///   1. Generic (pencil-beam): {dataDir}/machines/machine_{radiationMode}_{machineName}.json
///   2. Phase-space (folder):  {dataDir}/machines/{machineName}/{machineName}.json
///      (with .IAEAphsp/.IAEAheader files co-located in the same folder)
///
/// The loader auto-detects the schema from the JSON root keys.
class MachineLoader {
public:
    /// Load a machine by radiation mode and name.
    /// Tries folder-based path first, then flat file convention.
    /// @throws std::runtime_error if the file is not found or parsing fails
    static Machine load(const std::string& radiationMode,
                        const std::string& machineName,
                        const std::string& dataDir);

    /// Load a machine from an explicit JSON file path.
    /// @throws std::runtime_error if the file is not found or parsing fails
    static Machine loadFromFile(const std::string& filePath);

    /// Load a machine using the compiled-in default data directory.
    /// Equivalent to load(radiationMode, machineName, OPTIRAD_DATA_DIR).
    static Machine load(const std::string& radiationMode,
                        const std::string& machineName);

private:
    /// Load a machine from a "machine" JSON schema (phase-space / folder layout)
    static void loadPhaseSpaceMachine(const nlohmann::json& j,
                                       const std::string& jsonDir,
                                       Machine& machine);

    /// Load a machine from the traditional "meta"/"data"/"constraints" schema
    static void loadGenericMachine(const nlohmann::json& j,
                                    Machine& machine);

    /// Scan a directory for IAEA phase-space file pairs (.IAEAheader + .IAEAphsp)
    static std::vector<std::string> discoverPhaseSpaceFiles(const std::string& dir);
};

} // namespace optirad
