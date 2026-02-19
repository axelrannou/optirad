#pragma once

#include "core/Machine.hpp"
#include <string>

namespace optirad {

/// Loads machine data from JSON files.
///
/// Machine files follow the naming convention:
///   machine_{radiationMode}_{machineName}.json
///
/// and live under {dataDir}/machines/.
class MachineLoader {
public:
    /// Load a machine by radiation mode and name.
    /// Constructs the path: {dataDir}/machines/machine_{radiationMode}_{machineName}.json
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
};

} // namespace optirad
