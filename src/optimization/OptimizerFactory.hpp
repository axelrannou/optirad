#pragma once

#include "IOptimizer.hpp"
#include <memory>
#include <string>

namespace optirad {

class OptimizerFactory {
public:
    static std::unique_ptr<IOptimizer> create(const std::string& optimizerName);
};

} // namespace optirad
