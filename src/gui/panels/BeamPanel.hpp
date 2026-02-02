#pragma once

#include "IPanel.hpp"
#include "Plan.hpp"

namespace optirad {

class BeamPanel : public IPanel {
public:
    std::string getName() const override { return "Beams"; }
    void update() override;
    void render() override;

    void setPlan(Plan* plan);

private:
    Plan* m_plan = nullptr;
};

} // namespace optirad
