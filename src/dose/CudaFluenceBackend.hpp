#pragma once

#include <cstddef>
#include <memory>
#include <vector>

namespace optirad {

class FluenceDoseModel;
class Stf;
struct DoseObjectiveSpec;

class CudaFluenceBackend {
public:
    explicit CudaFluenceBackend(const FluenceDoseModel& model, size_t memoryBudgetBytes);
    ~CudaFluenceBackend();
    std::vector<double> computeDose(const std::vector<std::vector<float>>& fluence);
    std::vector<std::vector<float>> computeGradient(const std::vector<double>& voxelGradient);
    size_t memoryBytes() const;
    bool geometryResident() const;
    struct BixelLayout;
    std::shared_ptr<BixelLayout> createBixelLayout(const Stf& stf, double pixelMm, double boxMm);
    std::vector<double> computeBixelDose(BixelLayout& layout, const std::vector<double>& weights);
    void accumulateBixelGradient(BixelLayout& layout, const std::vector<double>& voxelGradient,
                                 std::vector<double>& gradient);
    bool evaluateObjectives(BixelLayout& layout, const std::vector<double>& weights,
        const std::vector<DoseObjectiveSpec>& objectives, double prescription, double threshold, double penalty,
        double& value, std::vector<double>& gradient, bool consistentHotspotGradient, bool needGradient);

private:
    struct Impl;
    std::unique_ptr<Impl> m_impl;
};

} // namespace optirad