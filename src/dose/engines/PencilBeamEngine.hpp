#pragma once

#include "../IDoseEngine.hpp"

namespace optirad {

class PencilBeamEngine : public IDoseEngine {
public:
    std::string getName() const override;

    DoseMatrix calculateDose(const Plan& plan, const Grid& grid) override;
    DoseInfluenceMatrix calculateDij(const Plan& plan, const Grid& grid) override;
};

} // namespace optirad
