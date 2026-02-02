#include "PencilBeamEngine.hpp"
#include "Logger.hpp"

namespace optirad {

std::string PencilBeamEngine::getName() const {
    return "PencilBeam";
}

DoseMatrix PencilBeamEngine::calculateDose(const Plan& plan, const Grid& grid) {
    Logger::info("Calculating dose with PencilBeam engine...");
    
    DoseMatrix dose;
    dose.setGrid(grid);
    dose.allocate();
    
    // TODO: Implement pencil beam dose calculation
    
    return dose;
}

DoseInfluenceMatrix PencilBeamEngine::calculateDij(const Plan& plan, const Grid& grid) {
    Logger::info("Calculating Dij matrix with PencilBeam engine...");
    
    DoseInfluenceMatrix dij;
    // TODO: Implement Dij calculation
    
    return dij;
}

} // namespace optirad
