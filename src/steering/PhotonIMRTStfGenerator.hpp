#pragma once

#include "IStfGenerator.hpp"
#include "../core/Beam.hpp"
#include "../core/Stf.hpp"
#include "../core/Machine.hpp"
#include "../geometry/MathUtils.hpp"
#include "../geometry/Grid.hpp"
#include "../geometry/Structure.hpp"
#include "../geometry/StructureSet.hpp"
#include "../geometry/VoxelDilation.hpp"
#include <memory>
#include <vector>
#include <set>
#include <algorithm>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

class PhotonIMRTStfGenerator : public IStfGenerator {
public:
    PhotonIMRTStfGenerator(double start = 0.0, double step = 60.0, double stop = 360.0,
                           double bixelWidth = 7.0, const std::array<double, 3>& iso = {0.0, 0.0, 0.0})
        : m_start(start), m_step(step), m_stop(stop), m_bixelWidth(bixelWidth), m_iso(iso) {}

    /// Set couch angle range (start/step/stop in degrees).
    /// Step <= 0 means a single couch angle equal to start for every beam.
    void setCouchAngles(double start, double step, double stop) {
        m_couchStart = start;
        m_couchStep = step;
        m_couchStop = stop;
        m_couchAnglesList.clear();
    }

    /// Set couch angles from an explicit paired list (like matRad).
    /// Takes precedence over start/step/stop if non-empty.
    void setCouchAngles(const std::vector<double>& angles) {
        m_couchAnglesList = angles;
    }

    /// Set gantry angles from an explicit list.
    /// Takes precedence over start/step/stop if non-empty.
    void setGantryAngles(const std::vector<double>& angles) {
        m_gantryAnglesList = angles;
    }

    /// Set the machine (used for SAD, SCD, energy)
    void setMachine(const Machine& machine) { m_machine = machine; m_hasMachine = true; }

    /// Set the field size [width, height] in mm (used only when no target voxels)
    void setFieldSize(const std::array<double, 2>& fieldSize) { m_fieldSize = fieldSize; }

    /// Set radiation mode label
    void setRadiationMode(const std::string& mode) { m_radiationMode = mode; }

    /// Set pre-computed target voxel world coordinates (LPS) for ray placement.
    /// These should already include any 3D margin expansion.
    void setTargetVoxelCoords(const std::vector<Vec3>& coords) { m_targetCoords = coords; }
    void setTargetVoxelCoords(std::vector<Vec3>&& coords) { m_targetCoords = std::move(coords); }

    /// Set CT resolution for bixel padding when bixelWidth < CT resolution
    void setCTResolution(const Vec3& res) { m_ctResolution = res; }

    /// Set the CT grid and structure set for 3D margin expansion.
    /// When both are provided, the generator performs proper 3D morphological
    /// dilation on the voxel grid before projecting onto the BEV plane.
    void setGrid(const Grid& grid) { m_grid = &grid; }
    void setStructureSet(const StructureSet& ss) { m_structureSet = &ss; }

    std::unique_ptr<StfProperties> generate() const override {
        auto props = std::make_unique<StfProperties>();
        props->bixelWidth = m_bixelWidth;
        if (!m_gantryAnglesList.empty()) {
            props->setGantryAngles(m_gantryAnglesList);
        } else {
            props->setGantryAngles(m_start, m_step, m_stop);
        }

        // Apply couch angles: explicit list takes precedence
        if (!m_couchAnglesList.empty()) {
            props->setCouchAngles(m_couchAnglesList);
        } else if (m_couchStep > 0.0) {
            props->setCouchAngles(m_couchStart, m_couchStep, m_couchStop);
        } else {
            props->setUniformCouchAngle(m_couchStart);
        }
        props->ensureConsistentAngles();

        props->setUniformIsoCenter(m_iso);
        return props;
    }

    /// Generate a full Stf object with beams and rays
    Stf generateStf() const {
        Stf stf;

        double sad = m_hasMachine ? m_machine.getSAD() : 1000.0;
        double scd = m_hasMachine ? m_machine.getSCD() : 500.0;
        double energy = m_hasMachine ? m_machine.getData().energy : 6.0;
        std::string machineName = m_hasMachine ? m_machine.getName() : "Generic";

        // Compute target world coords with 3D margin expansion if Grid/StructureSet available
        std::cout << "Computing target coordinates with 3D margin expansion..." << std::flush;
        std::vector<Vec3> targetCoords = computeTargetCoordsWithMargin();
        std::cout << " done (" << targetCoords.size() << " voxels)\n" << std::flush;
        
        bool hasTarget = !targetCoords.empty();
        Vec3 ctRes = m_ctResolution;

        // Collect gantry angles
        std::vector<double> gantryBase;
        if (!m_gantryAnglesList.empty()) {
            gantryBase = m_gantryAnglesList;
        } else {
            for (double angle = m_start; angle < m_stop; angle += m_step) {
                gantryBase.push_back(angle);
            }
        }

        // Collect couch angles and combine with gantry
        std::vector<double> angles;    // final gantry angles (after expansion)
        std::vector<double> couchAngles; // final couch angles (paired 1:1)

        if (!m_couchAnglesList.empty()) {
            // Explicit list: paired 1:1
            angles = gantryBase;
            couchAngles = m_couchAnglesList;
            // Resize if needed (pad with last value)
            if (couchAngles.size() != angles.size()) {
                double last = couchAngles.empty() ? 0.0 : couchAngles.back();
                couchAngles.resize(angles.size(), last);
            }
        } else if (m_couchStep > 0.0) {
            // Cartesian product: for each couch angle, ALL gantry angles (multi-arc)
            std::vector<double> couchRange;
            for (double a = m_couchStart; a < m_couchStop; a += m_couchStep) {
                couchRange.push_back(a);
            }
            if (couchRange.size() <= 1) {
                angles = gantryBase;
                couchAngles.assign(gantryBase.size(), couchRange.empty() ? m_couchStart : couchRange[0]);
            } else {
                for (double c : couchRange) {
                    for (double g : gantryBase) {
                        angles.push_back(g);
                        couchAngles.push_back(c);
                    }
                }
            }
        } else {
            // Single couch angle for all beams
            angles = gantryBase;
            couchAngles.assign(gantryBase.size(), m_couchStart);
        }

        // Generate beams in parallel
        std::vector<Beam> beams(angles.size());

#ifdef _OPENMP
        std::cout << "Generating " << angles.size() << " beams in parallel with OpenMP (" 
                  << omp_get_max_threads() << " threads)..." << std::flush;
        #pragma omp parallel for schedule(dynamic)
#else
        std::cout << "Generating " << angles.size() << " beams sequentially..." << std::flush;
#endif
        for (size_t i = 0; i < angles.size(); ++i) {
            double angle = angles[i];
            Beam& beam = beams[i];
            
            beam.setGantryAngle(angle);
            beam.setCouchAngle(couchAngles[i]);
            beam.setIsocenter(m_iso);
            beam.setBixelWidth(m_bixelWidth);
            beam.setRadiationMode(m_radiationMode);
            beam.setMachineName(machineName);
            beam.setSAD(sad);
            beam.setSCD(scd);

            if (hasTarget) {
                beam.generateRaysFromTarget(targetCoords, m_bixelWidth, ctRes);
            } else {
                beam.generateRays(m_bixelWidth, m_fieldSize);
            }

            beam.setAllRayEnergies(energy);
        }

        std::cout << " done\n" << std::flush;

        // Add beams to STF in order
        for (auto& beam : beams) {
            stf.addBeam(std::move(beam));
        }

        return stf;
    }

private:
    /// Compute target voxel world coordinates with 3D margin expansion.
    /// If Grid and StructureSet are available, performs proper 3D dilation.
    /// Otherwise falls back to the pre-set target coordinates.
    std::vector<Vec3> computeTargetCoordsWithMargin() const {
        // If no grid/structureSet, return pre-computed coords as-is
        if (!m_grid || !m_structureSet) {
            return m_targetCoords;
        }

        // Collect target and all-structure voxel indices
        std::set<size_t> targetVoxelSet;
        std::set<size_t> allVoxelSet;

        for (size_t i = 0; i < m_structureSet->getCount(); ++i) {
            const auto* structure = m_structureSet->getStructure(i);
            if (!structure) continue;

            const auto& voxels = structure->getVoxelIndices();
            allVoxelSet.insert(voxels.begin(), voxels.end());

            const auto& type = structure->getType();
            if (type == "PTV" || type == "GTV" || type == "CTV" || type == "TARGET") {
                targetVoxelSet.insert(voxels.begin(), voxels.end());
            }
        }

        if (targetVoxelSet.empty()) {
            return m_targetCoords;
        }

        // Compute margin: max(ctResolution, pbMargin) per axis
        // For photon IMRT: pbMargin = bixelWidth
        double pbMargin = m_bixelWidth;
        Vec3 spacing = m_grid->getSpacing();
        std::array<double, 3> margin = {
            std::max(spacing[0], pbMargin),
            std::max(spacing[1], pbMargin),
            std::max(spacing[2], pbMargin)
        };

        // Perform 3D morphological dilation
        std::vector<size_t> targetVec(targetVoxelSet.begin(), targetVoxelSet.end());
        auto dims = m_grid->getDimensions();
        std::vector<size_t> expandedIndices = dilateVoxels(
            targetVec, allVoxelSet, dims, {spacing[0], spacing[1], spacing[2]}, margin);

        // Convert expanded voxel indices to world coordinates
        size_t totalVoxels = dims[0] * dims[1] * dims[2];
        std::vector<Vec3> worldCoords;
        worldCoords.reserve(expandedIndices.size());

        for (size_t voxelIdx : expandedIndices) {
            if (voxelIdx >= totalVoxels) continue;

            // Column-major: index = row + col*ny + slice*ny*nx
            size_t row   = voxelIdx % dims[0];
            size_t col   = (voxelIdx / dims[0]) % dims[1];
            size_t slice = voxelIdx / (dims[0] * dims[1]);

            Vec3 ijk = {static_cast<double>(row), static_cast<double>(col), static_cast<double>(slice)};
            Vec3 worldPos = m_grid->voxelToPatient(ijk);
            worldCoords.push_back(worldPos);
        }

        return worldCoords;
    }

    double m_start;
    double m_step;
    double m_stop;
    double m_bixelWidth;
    std::array<double, 3> m_iso;
    std::array<double, 2> m_fieldSize = {100.0, 100.0};
    std::string m_radiationMode = "photons";
    Machine m_machine;
    bool m_hasMachine = false;
    std::vector<Vec3> m_targetCoords;  // pre-computed target world coordinates (fallback)
    Vec3 m_ctResolution = {1.0, 1.0, 1.0};
    const Grid* m_grid = nullptr;
    const StructureSet* m_structureSet = nullptr;

    // Couch angle configuration (paired 1:1 with gantry angles)
    double m_couchStart = 0.0;
    double m_couchStep = 0.0;   // 0 = single couch angle (m_couchStart) for all beams
    double m_couchStop = 0.0;
    std::vector<double> m_couchAnglesList;  // explicit list, takes precedence if non-empty
    std::vector<double> m_gantryAnglesList;  // explicit list, takes precedence if non-empty
};

} // namespace optirad
