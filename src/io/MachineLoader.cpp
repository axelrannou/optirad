#include "MachineLoader.hpp"
#include "utils/Logger.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad {

using json = nlohmann::json;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

Machine MachineLoader::load(const std::string& radiationMode,
                            const std::string& machineName,
                            const std::string& dataDir) {
    std::string filePath = dataDir + "/machines/machine_" + radiationMode + "_" + machineName + ".json";
    return loadFromFile(filePath);
}

Machine MachineLoader::load(const std::string& radiationMode,
                            const std::string& machineName) {
    return load(radiationMode, machineName, OPTIRAD_DATA_DIR);
}

Machine MachineLoader::loadFromFile(const std::string& filePath) {
    Logger::info("Loading machine from: " + filePath);

    std::ifstream file(filePath);
    if (!file.is_open()) {
        throw std::runtime_error("MachineLoader: cannot open file: " + filePath);
    }

    json j;
    try {
        j = json::parse(file);
    } catch (const json::parse_error& e) {
        throw std::runtime_error("MachineLoader: JSON parse error in " + filePath + ": " + e.what());
    }

    Machine machine;

    // ── Meta ──────────────────────────────────────────────────────────────
    if (j.contains("meta")) {
        const auto& jm = j["meta"];
        MachineMeta meta;
        if (jm.contains("radiationMode")) meta.radiationMode = jm["radiationMode"].get<std::string>();
        if (jm.contains("dataType"))      meta.dataType      = jm["dataType"].get<std::string>();
        if (jm.contains("created_on"))    meta.createdOn      = jm["created_on"].get<std::string>();
        if (jm.contains("created_by"))    meta.createdBy      = jm["created_by"].get<std::string>();
        if (jm.contains("description"))   meta.description    = jm["description"].get<std::string>();
        if (jm.contains("name"))          meta.name           = jm["name"].get<std::string>();
        if (jm.contains("SAD"))           meta.SAD            = jm["SAD"].get<double>();
        if (jm.contains("SCD"))           meta.SCD            = jm["SCD"].get<double>();
        machine.setMeta(meta);
    }

    // ── Data ─────────────────────────────────────────────────────────────
    if (j.contains("data")) {
        const auto& jd = j["data"];
        MachineData data;

        // Scalar fields
        if (jd.contains("energy"))           data.energy           = jd["energy"].get<double>();
        if (jd.contains("m"))                data.m                = jd["m"].get<double>();
        if (jd.contains("penumbraFWHMatIso")) data.penumbraFWHMatIso = jd["penumbraFWHMatIso"].get<double>();

        // Betas (array of 3)
        if (jd.contains("betas") && jd["betas"].is_array()) {
            const auto& jb = jd["betas"];
            for (size_t i = 0; i < std::min(jb.size(), size_t(3)); ++i) {
                data.betas[i] = jb[i].get<double>();
            }
        }

        // Primary fluence: array of [offAxis, value] pairs
        if (jd.contains("primaryFluence") && jd["primaryFluence"].is_array()) {
            data.primaryFluence.reserve(jd["primaryFluence"].size());
            for (const auto& entry : jd["primaryFluence"]) {
                if (entry.is_array() && entry.size() >= 2) {
                    data.primaryFluence.push_back({entry[0].get<double>(), entry[1].get<double>()});
                }
            }
        }

        // Kernel positions
        if (jd.contains("kernelPos") && jd["kernelPos"].is_array()) {
            data.kernelPos.reserve(jd["kernelPos"].size());
            for (const auto& v : jd["kernelPos"]) {
                data.kernelPos.push_back(v.get<double>());
            }
        }

        // Kernel entries: array of {SSD, kernel1, kernel2, kernel3}
        if (jd.contains("kernel") && jd["kernel"].is_array()) {
            data.kernel.reserve(jd["kernel"].size());
            for (const auto& jk : jd["kernel"]) {
                KernelEntry entry;
                if (jk.contains("SSD")) entry.SSD = jk["SSD"].get<double>();

                if (jk.contains("kernel1") && jk["kernel1"].is_array()) {
                    entry.kernel1.reserve(jk["kernel1"].size());
                    for (const auto& v : jk["kernel1"]) {
                        entry.kernel1.push_back(v.get<double>());
                    }
                }
                if (jk.contains("kernel2") && jk["kernel2"].is_array()) {
                    entry.kernel2.reserve(jk["kernel2"].size());
                    for (const auto& v : jk["kernel2"]) {
                        entry.kernel2.push_back(v.get<double>());
                    }
                }
                if (jk.contains("kernel3") && jk["kernel3"].is_array()) {
                    entry.kernel3.reserve(jk["kernel3"].size());
                    for (const auto& v : jk["kernel3"]) {
                        entry.kernel3.push_back(v.get<double>());
                    }
                }

                data.kernel.push_back(std::move(entry));
            }
        }

        machine.setData(data);

        Logger::info("Machine data loaded: " + std::to_string(data.primaryFluence.size()) +
                     " fluence entries, " + std::to_string(data.kernel.size()) +
                     " kernel entries, " + std::to_string(data.kernelPos.size()) +
                     " kernel positions");
    }

    // ── Constraints ──────────────────────────────────────────────────────
    if (j.contains("constraints")) {
        const auto& jc = j["constraints"];
        MachineConstraints constraints;

        if (jc.contains("gantryRotationSpeed") && jc["gantryRotationSpeed"].is_array()) {
            const auto& v = jc["gantryRotationSpeed"];
            if (v.size() >= 2) constraints.gantryRotationSpeed = {v[0].get<double>(), v[1].get<double>()};
        }
        if (jc.contains("leafSpeed") && jc["leafSpeed"].is_array()) {
            const auto& v = jc["leafSpeed"];
            if (v.size() >= 2) constraints.leafSpeed = {v[0].get<double>(), v[1].get<double>()};
        }
        if (jc.contains("monitorUnitRate") && jc["monitorUnitRate"].is_array()) {
            const auto& v = jc["monitorUnitRate"];
            if (v.size() >= 2) constraints.monitorUnitRate = {v[0].get<double>(), v[1].get<double>()};
        }

        machine.setConstraints(constraints);
    }

    Logger::info("Machine '" + machine.getName() + "' loaded successfully (" +
                 machine.getRadiationMode() + ", SAD=" + std::to_string(machine.getSAD()) + "mm)");
    return machine;
}

} // namespace optirad
