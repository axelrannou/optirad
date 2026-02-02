#pragma once

#include "IPanel.hpp"

namespace optirad {

class OptimizationPanel : public IPanel {
public:
    std::string getName() const override { return "Optimization"; }
    void update() override;
    void render() override;
};

} // namespace optirad
