#pragma once

#include "FluenceDoseModel.hpp"
#include "CudaFluenceBackend.hpp"
#include "IDoseOperator.hpp"
#include "core/Stf.hpp"

#include <memory>
#include <mutex>

namespace optirad {

/// Matrix-free IDoseOperator over a FluenceDoseModel: bixel weights are rasterized into fluence images per
/// beam, and the dose is evaluated by FFT convolution. Nothing proportional to the number of Dij entries is stored.
class FluenceDoseOperator : public IDoseOperator {
public:
    FluenceDoseOperator(std::shared_ptr<const FluenceDoseModel> model, std::shared_ptr<const Stf> stf);

    size_t getNumVoxels() const override { return m_model->numVoxels(); }
    size_t getNumBixels() const override { return m_offsets.back(); }
    std::vector<double> computeDose(const std::vector<double>& weights) const override;
    void accumulateTransposeProduct(const std::vector<double>& voxelGrad, std::vector<double>& grad) const override;
    double getMaxValue() const override;
    bool hasAllRows() const override { return m_model->computedRows().empty(); }
    bool isRowComputed(size_t voxel) const override {
        return m_model->computedRows().empty() || m_model->computedRows()[voxel] != 0;
    }
    DoseCalcOptions engineOptions() const override;
    bool evaluateObjectives(const std::vector<double>& weights, const std::vector<DoseObjectiveSpec>& objectives,
        double prescription, double threshold, double penalty, double& value, std::vector<double>& gradient,
        bool consistentHotspotGradient = false, bool needGradient = true) const override;

    const FluenceDoseModel& model() const { return *m_model; }

private:
    std::shared_ptr<const FluenceDoseModel> m_model;
    std::shared_ptr<const Stf> m_stf;
    std::vector<size_t> m_offsets;
    std::shared_ptr<CudaFluenceBackend::BixelLayout> m_cudaBixels;
    mutable std::once_flag m_maxOnce;
    mutable double m_maxValue = 0.0;
};

class ApertureDoseOperator : public IDoseOperator {
public:
    ApertureDoseOperator(std::shared_ptr<const IDoseOperator> bixelOperator, const Stf& stf,
                        std::vector<double> baseWeights);
    size_t getNumVoxels() const override { return m_bixelOperator->getNumVoxels(); }
    size_t getNumBixels() const override { return m_offsets.size() - 1; }
    std::vector<double> computeDose(const std::vector<double>& weights) const override;
    void accumulateTransposeProduct(const std::vector<double>& voxelGradient, std::vector<double>& gradient) const override;
    double getMaxValue() const override;
    bool hasAllRows() const override { return m_bixelOperator->hasAllRows(); }
    bool isRowComputed(size_t voxel) const override { return m_bixelOperator->isRowComputed(voxel); }
    DoseCalcOptions engineOptions() const override { return m_bixelOperator->engineOptions(); }
    bool evaluateObjectives(const std::vector<double>& weights, const std::vector<DoseObjectiveSpec>& objectives,
        double prescription, double threshold, double penalty, double& value, std::vector<double>& gradient,
        bool consistentHotspotGradient = false, bool needGradient = true) const override;
private:
    std::vector<double> expand(const std::vector<double>& weights) const;
    void collapse(const std::vector<double>& bixelGradient, std::vector<double>& gradient) const;
    std::shared_ptr<const IDoseOperator> m_bixelOperator;
    std::vector<double> m_baseWeights;
    std::vector<size_t> m_offsets;
};

} // namespace optirad
