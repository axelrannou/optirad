#pragma once

#include "IPanel.hpp"
#include "../AppState.hpp"

namespace optirad {

class DeliveryPanel : public IPanel {
public:
    explicit DeliveryPanel(GuiAppState& state) : m_state(state) {}
    std::string getName() const override { return "Delivery"; }
    void update() override {}
    void render() override;
    static void renderControls(DeliveryPlayback& delivery);

private:
    GuiAppState& m_state;
};

} // namespace optirad
