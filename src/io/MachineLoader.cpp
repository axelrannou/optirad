#include "MachineLoader.hpp"
#include "utils/Logger.hpp"

#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>
#include <filesystem>
#include <algorithm>

#ifndef OPTIRAD_DATA_DIR
#define OPTIRAD_DATA_DIR "."
#endif

namespace optirad {

using json = nlohmann::json;
namespace fs = std::filesystem;

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

Machine MachineLoader::load(const std::string& radiationMode,
                            const std::string& machineName,
                            const std::string& dataDir) {
    // Try folder-based path first: {dataDir}/machines/{machineName}/{machineName}.json
    std::string folderPath = dataDir + "/machines/" + machineName + "/" + machineName + ".json";
    if (fs::exists(folderPath)) {
        Logger::info("MachineLoader: found folder-based machine at " + folderPath);
        return loadFromFile(folderPath);
    }

    // Fallback to flat file: {dataDir}/machines/machine_{radiationMode}_{machineName}.json
    std::string flatPath = dataDir + "/machines/machine_" + radiationMode + "_" + machineName + ".json";
    return loadFromFile(flatPath);
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

    // Auto-detect schema: "machine" key = phase-space, "meta" key = generic
    if (j.contains("machine")) {
        // Get the directory containing the JSON file for IAEA file discovery
        std::string jsonDir = fs::path(filePath).parent_path().string();
        loadPhaseSpaceMachine(j, jsonDir, machine);
    } else {
        loadGenericMachine(j, machine);
    }

    Logger::info("Machine '" + machine.getName() + "' loaded successfully (" +
                 machine.getRadiationMode() + ", SAD=" + std::to_string(machine.getSAD()) +
                 "mm, type=" + (machine.isPhaseSpace() ? "PhaseSpace" : "Generic") + ")");
    return machine;
}

// ---------------------------------------------------------------------------
// Generic machine loader (traditional "meta"/"data"/"constraints" schema)
// ---------------------------------------------------------------------------

void MachineLoader::loadGenericMachine(const json& j, Machine& machine) {
    // ── Meta ──────────────────────────────────────────────────────────────
    if (j.contains("meta")) {
        const auto& jm = j["meta"];
        MachineMeta meta;
        meta.machineType = MachineType::Generic;
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

        if (jd.contains("energy"))           data.energy           = jd["energy"].get<double>();
        if (jd.contains("m"))                data.m                = jd["m"].get<double>();
        if (jd.contains("penumbraFWHMatIso")) data.penumbraFWHMatIso = jd["penumbraFWHMatIso"].get<double>();

        if (jd.contains("betas") && jd["betas"].is_array()) {
            const auto& jb = jd["betas"];
            for (size_t i = 0; i < std::min(jb.size(), size_t(3)); ++i) {
                data.betas[i] = jb[i].get<double>();
            }
        }

        if (jd.contains("primaryFluence") && jd["primaryFluence"].is_array()) {
            data.primaryFluence.reserve(jd["primaryFluence"].size());
            for (const auto& entry : jd["primaryFluence"]) {
                if (entry.is_array() && entry.size() >= 2) {
                    data.primaryFluence.push_back({entry[0].get<double>(), entry[1].get<double>()});
                }
            }
        }

        if (jd.contains("kernelPos") && jd["kernelPos"].is_array()) {
            data.kernelPos.reserve(jd["kernelPos"].size());
            for (const auto& v : jd["kernelPos"]) {
                data.kernelPos.push_back(v.get<double>());
            }
        }

        if (jd.contains("kernel") && jd["kernel"].is_array()) {
            data.kernel.reserve(jd["kernel"].size());
            for (const auto& jk : jd["kernel"]) {
                KernelEntry entry;
                if (jk.contains("SSD")) entry.SSD = jk["SSD"].get<double>();

                if (jk.contains("kernel1") && jk["kernel1"].is_array()) {
                    entry.kernel1.reserve(jk["kernel1"].size());
                    for (const auto& v : jk["kernel1"]) entry.kernel1.push_back(v.get<double>());
                }
                if (jk.contains("kernel2") && jk["kernel2"].is_array()) {
                    entry.kernel2.reserve(jk["kernel2"].size());
                    for (const auto& v : jk["kernel2"]) entry.kernel2.push_back(v.get<double>());
                }
                if (jk.contains("kernel3") && jk["kernel3"].is_array()) {
                    entry.kernel3.reserve(jk["kernel3"].size());
                    for (const auto& v : jk["kernel3"]) entry.kernel3.push_back(v.get<double>());
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
}

// ---------------------------------------------------------------------------
// Phase-space machine loader (folder-based "machine" schema)
// ---------------------------------------------------------------------------

void MachineLoader::loadPhaseSpaceMachine(const json& j,
                                           const std::string& jsonDir,
                                           Machine& machine) {
    const auto& jm = j["machine"];

    // ── Meta ──
    MachineMeta meta;
    meta.machineType = MachineType::PhaseSpace;
    meta.radiationMode = "photons";

    if (jm.contains("name")) meta.name = jm["name"].get<std::string>();
    if (jm.contains("type")) {
        std::string type = jm["type"].get<std::string>();
        meta.radiationMode = (type == "photon" || type == "photons") ? "photons" : type;
    }
    if (jm.contains("source_to_isocenter_distance_cm")) {
        meta.SAD = jm["source_to_isocenter_distance_cm"].get<double>() * 10.0; // cm → mm
    }

    machine.setMeta(meta);

    // ── Geometry ──
    MachineGeometry geom;

    // Jaws (cm → mm)
    if (jm.contains("jaws")) {
        const auto& jj = jm["jaws"];
        if (jj.contains("X1_min_cm")) geom.jawX1Min = jj["X1_min_cm"].get<double>() * 10.0;
        if (jj.contains("X1_max_cm")) geom.jawX1Max = jj["X1_max_cm"].get<double>() * 10.0;
        if (jj.contains("X2_min_cm")) geom.jawX2Min = jj["X2_min_cm"].get<double>() * 10.0;
        if (jj.contains("X2_max_cm")) geom.jawX2Max = jj["X2_max_cm"].get<double>() * 10.0;
        if (jj.contains("Y1_min_cm")) geom.jawY1Min = jj["Y1_min_cm"].get<double>() * 10.0;
        if (jj.contains("Y1_max_cm")) geom.jawY1Max = jj["Y1_max_cm"].get<double>() * 10.0;
        if (jj.contains("Y2_min_cm")) geom.jawY2Min = jj["Y2_min_cm"].get<double>() * 10.0;
        if (jj.contains("Y2_max_cm")) geom.jawY2Max = jj["Y2_max_cm"].get<double>() * 10.0;
        if (jj.contains("default_field_size_cm") && jj["default_field_size_cm"].is_array()) {
            const auto& fs = jj["default_field_size_cm"];
            if (fs.size() >= 2) {
                geom.defaultFieldSize = {fs[0].get<double>() * 10.0, fs[1].get<double>() * 10.0};
            }
        }
    }

    // Collimator
    if (jm.contains("collimator")) {
        const auto& jc = jm["collimator"];
        if (jc.contains("min_deg") && jc.contains("max_deg")) {
            geom.collimatorRange = {jc["min_deg"].get<double>(), jc["max_deg"].get<double>()};
        }
        if (jc.contains("default_deg")) {
            geom.defaultCollimatorAngle = jc["default_deg"].get<double>();
        }
    }

    // Couch
    if (jm.contains("couch")) {
        const auto& jc = jm["couch"];
        if (jc.contains("min_deg") && jc.contains("max_deg")) {
            geom.couchRange = {jc["min_deg"].get<double>(), jc["max_deg"].get<double>()};
        }
        if (jc.contains("default_deg")) {
            geom.defaultCouchAngle = jc["default_deg"].get<double>();
        }
    }

    // MLC
    if (jm.contains("MLC")) {
        const auto& jmlc = jm["MLC"];
        if (jmlc.contains("type")) geom.mlcType = jmlc["type"].get<std::string>();
        if (jmlc.contains("num_leaves")) geom.numLeaves = jmlc["num_leaves"].get<int>();
        if (jmlc.contains("max_travel_mm")) geom.maxLeafTravel = jmlc["max_travel_mm"].get<double>();
        if (jmlc.contains("interdigitation")) geom.interdigitation = jmlc["interdigitation"].get<bool>();
        if (jmlc.contains("leaf_width_mm") && jmlc["leaf_width_mm"].is_array()) {
            for (const auto& w : jmlc["leaf_width_mm"]) {
                geom.leafWidths.push_back(w.get<double>());
            }
        }
    }

    // Beam energy / dose rate
    if (jm.contains("beam_energy_MV")) geom.beamEnergyMV = jm["beam_energy_MV"].get<double>();
    if (jm.contains("dose_rate_MU_per_min")) geom.doseRateMUPerMin = jm["dose_rate_MU_per_min"].get<double>();
    if (jm.contains("aperture_sampling_resolution_mm")) {
        geom.apertureSamplingResolutionMm = jm["aperture_sampling_resolution_mm"].get<double>();
    }
    if (jm.contains("max_field_size_cm") && jm["max_field_size_cm"].is_array()) {
        const auto& mfs = jm["max_field_size_cm"];
        if (mfs.size() >= 2) {
            geom.defaultFieldSize = {mfs[0].get<double>() * 10.0, mfs[1].get<double>() * 10.0};
        }
    }

    // Discover IAEA phase-space files in the same directory
    geom.phaseSpaceDir = jsonDir;
    geom.phaseSpaceFileNames = discoverPhaseSpaceFiles(jsonDir);
    geom.numPhaseSpaceFiles = static_cast<int>(geom.phaseSpaceFileNames.size());

    machine.setGeometry(geom);

    Logger::info("Phase-space machine: " + meta.name +
                 ", SAD=" + std::to_string(meta.SAD) + "mm" +
                 ", MLC=" + geom.mlcType +
                 ", " + std::to_string(geom.numPhaseSpaceFiles) + " PSF file(s) found");
}

// ---------------------------------------------------------------------------
// IAEA file discovery
// ---------------------------------------------------------------------------

std::vector<std::string> MachineLoader::discoverPhaseSpaceFiles(const std::string& dir) {
    std::vector<std::string> baseNames;

    if (!fs::exists(dir) || !fs::is_directory(dir)) {
        Logger::warn("MachineLoader: Phase-space directory does not exist: " + dir);
        return baseNames;
    }

    for (const auto& entry : fs::directory_iterator(dir)) {
        if (!entry.is_regular_file()) continue;
        std::string filename = entry.path().filename().string();
        // Look for .IAEAheader files
        if (filename.size() > 11 && filename.substr(filename.size() - 11) == ".IAEAheader") {
            std::string baseName = filename.substr(0, filename.size() - 11);
            // Verify the companion .IAEAphsp file exists
            std::string phspPath = dir + "/" + baseName + ".IAEAphsp";
            if (fs::exists(phspPath)) {
                baseNames.push_back(baseName);
            } else {
                Logger::warn("MachineLoader: header without companion phsp: " + filename);
            }
        }
    }

    // Sort for deterministic ordering
    std::sort(baseNames.begin(), baseNames.end());

    return baseNames;
}

} // namespace optirad
