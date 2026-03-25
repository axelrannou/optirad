#include "PhaseSpaceBuilder.hpp"
#include "utils/Logger.hpp"
#include <stdexcept>

namespace optirad {

std::vector<std::shared_ptr<PhaseSpaceBeamSource>> PhaseSpaceBuilder::build(
    const Plan& plan,
    const PhaseSpaceBuildOptions& options) {

    if (!plan.getMachine().isPhaseSpace()) {
        throw std::runtime_error("PhaseSpaceBuilder: machine is not a phase-space machine");
    }

    const auto& stfProps = plan.getStfProperties();
    const auto& gantryAngles = stfProps.gantryAngles;
    const auto& couchAngles = stfProps.couchAngles;
    std::array<double, 3> iso = {0.0, 0.0, 0.0};
    if (!stfProps.isoCenters.empty()) {
        iso = stfProps.isoCenters[0];
    }

    const int numBeams = static_cast<int>(gantryAngles.size());
    std::vector<std::shared_ptr<PhaseSpaceBeamSource>> sources(numBeams);

    Logger::info("PhaseSpaceBuilder: building " + std::to_string(numBeams) +
                 " beam sources (maxParticles=" + std::to_string(options.maxParticles) + ")");

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < numBeams; ++i) {
        auto source = std::make_shared<PhaseSpaceBeamSource>();
        double couch = (static_cast<size_t>(i) < couchAngles.size())
                           ? couchAngles[i] : 0.0;
        source->configure(plan.getMachine(), gantryAngles[i],
                          options.collimatorAngle, couch, iso);
        source->build(options.maxParticles, options.vizSampleSize);
        sources[i] = std::move(source);
    }

    Logger::info("PhaseSpaceBuilder: built " + std::to_string(numBeams) + " beam sources");
    return sources;
}

} // namespace optirad
