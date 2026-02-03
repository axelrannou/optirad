#pragma once

#include "IPanel.hpp"
#include <vector>
#include <string>

namespace optirad {

class LogPanel : public IPanel {
public:
    std::string getName() const override { return "Log"; }
    void update() override;
    void render() override;

    void addLog(const std::string& message);
    void clear();

private:
    std::vector<std::string> m_logs;
};

} // namespace optirad
