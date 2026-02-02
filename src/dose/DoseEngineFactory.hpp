#pragma once

#include "IDoseEngine.hpp"
#include <memory>
#include <string>

namespace optirad {

class DoseEngineFactory {
public:
    static std::unique_ptr<IDoseEngine> create(const std::string& engineName);
};

} // namespace optirad
