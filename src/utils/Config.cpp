#include "Config.hpp"

namespace optirad {

bool Config::load(const std::string& filepath) {
    // TODO: Load JSON config
    return true;
}

bool Config::save(const std::string& filepath) const {
    // TODO: Save JSON config
    return true;
}

std::string Config::getDoseEngine() const { return m_doseEngine; }
std::string Config::getOptimizer() const { return m_optimizer; }

} // namespace optirad
