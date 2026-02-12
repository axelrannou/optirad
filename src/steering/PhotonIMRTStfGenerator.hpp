#pragma once

#include "IStfGenerator.hpp"
#include <memory>

namespace optirad {

class PhotonIMRTStfGenerator : public IStfGenerator {
public:
    PhotonIMRTStfGenerator(double start = 0.0, double step = 60.0, double stop = 360.0, double bixelWidth = 7.0, const std::array<double, 3>& iso = {0.0, 0.0, 0.0})
        : m_start(start), m_step(step), m_stop(stop), m_bixelWidth(bixelWidth), m_iso(iso) {}

    std::unique_ptr<StfProperties> generate() const override {
        auto props = std::make_unique<StfProperties>();
        props->bixelWidth = m_bixelWidth;
        props->setGantryAngles(m_start, m_step, m_stop);
        props->setUniformIsoCenter(m_iso);
        return props;
    }

private:
    double m_start;
    double m_step;
    double m_stop;
    double m_bixelWidth;
    std::array<double, 3> m_iso;
};

} // namespace optirad
