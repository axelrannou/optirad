#pragma once

#include "StfProperties.hpp"
#include <memory>

namespace optirad {

class IStfGenerator {
public:
    virtual ~IStfGenerator() = default;
    virtual std::unique_ptr<StfProperties> generate() const = 0;
};

} // namespace optirad
