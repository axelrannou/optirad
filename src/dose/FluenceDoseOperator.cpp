#include "FluenceDoseOperator.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace optirad {

FluenceDoseOperator::FluenceDoseOperator(std::shared_ptr<const FluenceDoseModel> model,
                                         std::shared_ptr<const Stf> stf)
    : m_model(std::move(model)), m_stf(std::move(stf)) {
    if (!m_model || !m_stf || m_model->numBeams() != m_stf->getCount())
        throw std::invalid_argument("FluenceDoseOperator: model and STF must have the same beams");
    m_offsets.assign(m_stf->getCount() + 1, 0);
    for (size_t b = 0; b < m_stf->getCount(); ++b)
        m_offsets[b + 1] = m_offsets[b] + m_stf->getBeam(b)->getTotalNumOfBixels();
    if (m_model->hasPrimaryFluence()) {
        m_bixelScale.reserve(m_offsets.back());
        for (size_t b = 0; b < m_stf->getCount(); ++b) {
            const Beam& beam = *m_stf->getBeam(b);
            for (size_t ri = 0; ri < beam.getNumOfRays(); ++ri) {
                const auto* ray = beam.getRay(ri);
                if (!ray) continue;
                const auto p = ray->getRayPosBev();
                m_bixelScale.insert(m_bixelScale.end(), ray->getNumOfBixels(),
                                    m_model->relativeFluence(std::hypot(p[0], p[2])));
            }
        }
        if (m_bixelScale.size() != m_offsets.back()) m_bixelScale.clear();
    }
#ifdef OPTIRAD_HAS_CUDA
    if (m_model->m_cuda) m_cudaBixels = m_model->m_cuda->createBixelLayout(*m_stf, m_model->m_pixel, m_model->m_boxMm);
#endif
}

std::vector<double> FluenceDoseOperator::scaled(const std::vector<double>& weights) const {
    std::vector<double> out(weights);
    for (size_t i = 0; i < out.size() && i < m_bixelScale.size(); ++i) out[i] *= m_bixelScale[i];
    return out;
}

void FluenceDoseOperator::scaleGradient(std::vector<double>& grad, size_t begin) const {
    for (size_t i = begin; i < grad.size() && i < m_bixelScale.size(); ++i) grad[i] *= m_bixelScale[i];
}

std::vector<double> FluenceDoseOperator::computeDose(const std::vector<double>& rawWeights) const {
    if (rawWeights.size() != getNumBixels()) throw std::invalid_argument("FluenceDoseOperator: bad weight count");
    const std::vector<double> scaledWeights = m_bixelScale.empty() ? std::vector<double>{} : scaled(rawWeights);
    const std::vector<double>& weights = m_bixelScale.empty() ? rawWeights : scaledWeights;
#ifdef OPTIRAD_HAS_CUDA
    if (m_cudaBixels) return m_model->m_cuda->computeBixelDose(*m_cudaBixels, weights);
#endif
    const size_t nb = m_stf->getCount();
    std::vector<std::vector<float>> images(nb);
#pragma omp parallel for schedule(dynamic, 1)
    for (long b = 0; b < static_cast<long>(nb); ++b) {
        const size_t bi = static_cast<size_t>(b);
        const double* w = weights.data() + m_offsets[bi];
        bool any = false;
        for (size_t i = 0; i < m_offsets[bi + 1] - m_offsets[bi] && !any; ++i) any = w[i] != 0.0;
        if (any) m_model->rasterizeBixels(*m_stf->getBeam(bi), w, images[bi]);
    }
    return m_model->computeDose(images);
}

void FluenceDoseOperator::accumulateTransposeProduct(const std::vector<double>& voxelGrad,
                                                     std::vector<double>& grad) const {
    if (grad.size() != getNumBixels()) throw std::invalid_argument("FluenceDoseOperator: bad gradient size");
    if (voxelGrad.size() != getNumVoxels()) throw std::invalid_argument("FluenceDoseOperator: bad voxel gradient size");
    // d(dose)/d(weight) carries the bixel's fluence factor: gather into a zeroed buffer, scale, then add.
    std::vector<double> local;
    std::vector<double>& target = m_bixelScale.empty() ? grad : (local.assign(grad.size(), 0.0), local);
#ifdef OPTIRAD_HAS_CUDA
    if (m_cudaBixels) {
        m_model->m_cuda->accumulateBixelGradient(*m_cudaBixels, voxelGrad, target);
    } else
#endif
    {
        const auto images = m_model->computeGradient(voxelGrad);
#pragma omp parallel for schedule(dynamic, 1)
        for (long b = 0; b < static_cast<long>(images.size()); ++b) {
            const size_t bi = static_cast<size_t>(b);
            m_model->gatherBixels(*m_stf->getBeam(bi), images[bi], target.data() + m_offsets[bi]);
        }
    }
    if (!m_bixelScale.empty())
        for (size_t i = 0; i < grad.size(); ++i) grad[i] += local[i] * m_bixelScale[i];
}

double FluenceDoseOperator::getMaxValue() const {
    std::call_once(m_maxOnce, [this] {
        // The bixel nearest the axis of a few evenly spaced beams is the strongest single bixel.
        const size_t nb = m_stf->getCount();
        const size_t samples = std::min<size_t>(nb, 8);
        std::vector<double> dose(m_model->numVoxels());
        for (size_t s = 0; s < samples; ++s) {
            const size_t b = s * nb / samples;
            const Beam& beam = *m_stf->getBeam(b);
            std::vector<double> w(beam.getTotalNumOfBixels(), 0.0);
            size_t best = 0, bixel = 0;
            double bestR = 1e30;
            for (size_t ri = 0; ri < beam.getNumOfRays(); ++ri) {
                const auto p = beam.getRay(ri)->getRayPosBev();
                const double r = std::hypot(p[0], p[2]);
                if (r < bestR) { bestR = r; best = bixel; }
                bixel += beam.getRay(ri)->getNumOfBixels();
            }
            const double probeWeight = m_bixelScale.empty() ? 1.0 : m_bixelScale[m_offsets[b] + best];
            w[best] = probeWeight;
#ifdef OPTIRAD_HAS_CUDA
            if (m_cudaBixels) {
                std::vector<double> probe(getNumBixels(), 0.0);
                probe[m_offsets[b] + best] = 1.0;  // computeDose applies the fluence factor itself
                const auto values = computeDose(probe);
                m_maxValue = std::max(m_maxValue, *std::max_element(values.begin(), values.end()));
                continue;
            }
#endif
            std::vector<float> image;
            m_model->rasterizeBixels(beam, w.data(), image);
            std::fill(dose.begin(), dose.end(), 0.0);
            m_model->addBeamDose(b, image, dose.data());
            m_maxValue = std::max(m_maxValue, *std::max_element(dose.begin(), dose.end()));
        }
    });
    return m_maxValue;
}

DoseCalcOptions FluenceDoseOperator::engineOptions() const {
    DoseCalcOptions o;
    o.divergentBixels = true;
    o.gpuRayTracing = m_model->usesCuda();
    return o;
}

bool FluenceDoseOperator::evaluateObjectives(const std::vector<double>& weights,
    const std::vector<DoseObjectiveSpec>& objectives, double& value, std::vector<double>& gradient,
    bool needGradient) const {
#ifdef OPTIRAD_HAS_CUDA
    if (m_cudaBixels) {
        if (m_bixelScale.empty())
            return m_model->m_cuda->evaluateObjectives(*m_cudaBixels, weights, objectives, value, gradient, needGradient);
        const bool ok = m_model->m_cuda->evaluateObjectives(*m_cudaBixels, scaled(weights), objectives, value,
            gradient, needGradient);
        if (ok && needGradient) scaleGradient(gradient, 0);
        return ok;
    }
#else
    (void)weights; (void)objectives; (void)value; (void)gradient; (void)needGradient;
#endif
    return false;
}

ApertureDoseOperator::ApertureDoseOperator(std::shared_ptr<const IDoseOperator> bixelOperator,
    const Stf& stf, std::vector<double> baseWeights)
    : m_bixelOperator(std::move(bixelOperator)), m_baseWeights(std::move(baseWeights)), m_offsets(1, 0) {
    if (!m_bixelOperator || m_baseWeights.size() != m_bixelOperator->getNumBixels() ||
        m_baseWeights.size() != stf.getTotalNumOfBixels())
        throw std::invalid_argument("ApertureDoseOperator: incompatible bixel dimensions");
    for (const auto& beam : stf.getBeams()) m_offsets.push_back(m_offsets.back() + beam.getTotalNumOfBixels());
}

std::vector<double> ApertureDoseOperator::expand(const std::vector<double>& weights) const {
    if (weights.size() != getNumBixels()) throw std::invalid_argument("ApertureDoseOperator: bad aperture weight count");
    auto expanded = m_baseWeights;
    for (size_t aperture = 0; aperture < weights.size(); ++aperture)
        for (size_t bixel = m_offsets[aperture]; bixel < m_offsets[aperture + 1]; ++bixel)
            expanded[bixel] *= weights[aperture];
    return expanded;
}

void ApertureDoseOperator::collapse(const std::vector<double>& bixelGradient, std::vector<double>& gradient) const {
    if (gradient.size() != getNumBixels()) throw std::invalid_argument("ApertureDoseOperator: bad aperture gradient count");
    for (size_t aperture = 0; aperture < gradient.size(); ++aperture)
        for (size_t bixel = m_offsets[aperture]; bixel < m_offsets[aperture + 1]; ++bixel)
            gradient[aperture] += m_baseWeights[bixel] * bixelGradient[bixel];
}

std::vector<double> ApertureDoseOperator::computeDose(const std::vector<double>& weights) const {
    return m_bixelOperator->computeDose(expand(weights));
}

void ApertureDoseOperator::accumulateTransposeProduct(const std::vector<double>& voxelGradient,
    std::vector<double>& gradient) const {
    std::vector<double> bixelGradient(m_baseWeights.size(), 0.0);
    m_bixelOperator->accumulateTransposeProduct(voxelGradient, bixelGradient);
    collapse(bixelGradient, gradient);
}

double ApertureDoseOperator::getMaxValue() const {
    double maximum = 0.0;
    std::vector<double> weights(getNumBixels(), 0.0);
    for (size_t aperture = 0; aperture < weights.size(); ++aperture) {
        weights[aperture] = 1.0;
        const auto dose = computeDose(weights);
        if (!dose.empty()) maximum = std::max(maximum, *std::max_element(dose.begin(), dose.end()));
        weights[aperture] = 0.0;
    }
    return maximum;
}

bool ApertureDoseOperator::evaluateObjectives(const std::vector<double>& weights,
    const std::vector<DoseObjectiveSpec>& objectives, double& value, std::vector<double>& gradient,
    bool needGradient) const {
    std::vector<double> bixelGradient;
    if (!m_bixelOperator->evaluateObjectives(expand(weights), objectives, value, bixelGradient, needGradient))
        return false;
    gradient.assign(needGradient ? getNumBixels() : 0, 0.0);
    if (needGradient) collapse(bixelGradient, gradient);
    return true;
}

} // namespace optirad
