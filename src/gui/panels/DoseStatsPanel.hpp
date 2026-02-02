#pragma once

#include "IPanel.hpp"

namespace optirad {

class DoseStatsPanel : public IPanel {
public:
    std::string getName() const override { return "Dose Statistics"; }
    void update() override;
    void render() override;
};

} // namespace optirad
