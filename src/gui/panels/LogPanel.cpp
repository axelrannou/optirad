#include "LogPanel.hpp"

namespace optirad {

void LogPanel::update() {}

void LogPanel::render() {
    if (!m_visible) return;
    
    // TODO: ImGui::Begin("Log");
    // for (const auto& log : m_logs) {
    //     ImGui::TextUnformatted(log.c_str());
    // }
    // if (ImGui::Button("Clear")) clear();
    // ImGui::End();
}

void LogPanel::addLog(const std::string& message) { m_logs.push_back(message); }
void LogPanel::clear() { m_logs.clear(); }

} // namespace optirad
