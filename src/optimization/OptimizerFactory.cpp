#include "OptimizerFactory.hpp"
#include "optimizers/LBFGSOptimizer.hpp"
#include <stdexcept>

namespace optirad {

std::unique_ptr<IOptimizer> OptimizerFactory::create(const std::string& optimizerName) {
    if (optimizerName == "LBFGS") {
        return std::make_unique<LBFGSOptimizer>();
    }
    // Add more optimizers here:
    // else if (optimizerName == "IPOPT") {
    //     return std::make_unique<IPOPTOptimizer>();
    // }
    
    throw std::runtime_error("Unknown optimizer: " + optimizerName);
}

} // namespace optirad
