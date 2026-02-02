#pragma once

#include <string>

namespace optirad {

class Config {
public:
    bool load(const std::string& filepath);
    bool save(const std::string& filepath) const;

    std::string getDoseEngine() const;
    std::string getOptimizer() const;

private:
    std::string m_doseEngine = "PencilBeam";
    std::string m_optimizer = "LBFGS";
};

} // namespace optirad
