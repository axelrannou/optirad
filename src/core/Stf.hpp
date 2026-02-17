#pragma once

#include "../core/Beam.hpp"
#include "../core/Ray.hpp"
#include <vector>
#include <memory>
#include <string>

namespace optirad {

/**
 * Stf (Steering File) - manages the collection of beams for a treatment plan.
 * Analogous to a 1×N struct array in external tools, where each element is a beam
 * containing angular parameters, geometry, and a list of rays.
 */
class Stf {
public:
    Stf() = default;

    // Beam management
    void addBeam(const Beam& beam);
    void addBeam(Beam&& beam);
    
    const Beam* getBeam(size_t index) const;
    Beam* getBeam(size_t index);
    
    size_t getCount() const { return m_beams.size(); }
    bool isEmpty() const { return m_beams.empty(); }
    
    const std::vector<Beam>& getBeams() const { return m_beams; }
    std::vector<Beam>& beams() { return m_beams; }
    
    // Clear all beams
    void clear() { m_beams.clear(); }
    
    // Computed properties across all beams
    size_t getTotalNumOfRays() const;
    size_t getTotalNumOfBixels() const;
    
    // Utility methods
    void computeAllSourcePoints();
    void generateAllRays(double bixelWidth, const std::array<double, 2>& fieldSize);
    
    // Summary
    void printSummary() const;

private:
    std::vector<Beam> m_beams;
};

} // namespace optirad
