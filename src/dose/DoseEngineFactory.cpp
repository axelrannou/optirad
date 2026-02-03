#include "DoseEngineFactory.hpp"
#include "engines/PencilBeamEngine.hpp"
#include <stdexcept>

namespace optirad {

std::unique_ptr<IDoseEngine> DoseEngineFactory::create(const std::string& engineName) {
    if (engineName == "PencilBeam") {
        return std::make_unique<PencilBeamEngine>();
    }
    // Add more engines here:
    // else if (engineName == "MonteCarlo") {
    //     return std::make_unique<MonteCarloEngine>();
    // }
    
    throw std::runtime_error("Unknown dose engine: " + engineName);
}

} // namespace optirad
