#include "FluenceDoseModel.hpp"
#include "CudaFluenceBackend.hpp"
#include "utils/Logger.hpp"

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace optirad {

namespace {

struct C { double r, i; };

class Fft2 {
public:
    explicit Fft2(size_t n) : m_n(n), m_rev(n), m_tw(n / 2) {
        size_t bits = 0;
        while ((size_t(1) << bits) < n) ++bits;
        for (size_t i = 0; i < n; ++i) {
            size_t r = 0;
            for (size_t b = 0; b < bits; ++b) if (i & (size_t(1) << b)) r |= size_t(1) << (bits - 1 - b);
            m_rev[i] = r;
        }
        for (size_t k = 0; k < n / 2; ++k) {
            const double a = -2.0 * M_PI * static_cast<double>(k) / static_cast<double>(n);
            m_tw[k] = {std::cos(a), std::sin(a)};
        }
    }

    void forward(C* d) const { rows(d); transpose(d); rows(d); transpose(d); }

    void inverse(C* d) const {
        const size_t total = m_n * m_n;
        for (size_t i = 0; i < total; ++i) d[i].i = -d[i].i;
        forward(d);
        const double s = 1.0 / static_cast<double>(total);
        for (size_t i = 0; i < total; ++i) { d[i].r *= s; d[i].i *= -s; }
    }

private:
    void rows(C* d) const { for (size_t r = 0; r < m_n; ++r) fft1(d + r * m_n); }

    void transpose(C* d) const {
        for (size_t r = 0; r < m_n; ++r)
            for (size_t c = r + 1; c < m_n; ++c) std::swap(d[r * m_n + c], d[c * m_n + r]);
    }

    void fft1(C* a) const {
        const size_t n = m_n;
        for (size_t i = 0; i < n; ++i) if (i < m_rev[i]) std::swap(a[i], a[m_rev[i]]);
        for (size_t len = 2; len <= n; len <<= 1) {
            const size_t half = len / 2, step = n / len;
            for (size_t i = 0; i < n; i += len) {
                for (size_t j = 0; j < half; ++j) {
                    const C w = m_tw[j * step];
                    const C v = a[i + j + half];
                    const C t = {v.r * w.r - v.i * w.i, v.r * w.i + v.i * w.r};
                    const C u = a[i + j];
                    a[i + j] = {u.r + t.r, u.i + t.i};
                    a[i + j + half] = {u.r - t.r, u.i - t.i};
                }
            }
        }
    }

    size_t m_n;
    std::vector<size_t> m_rev;
    std::vector<C> m_tw;
};

struct Workspace {
    size_t n = 0;
    std::vector<C> a, b;
    std::vector<float> g[3];
    std::unique_ptr<Fft2> fft;

    void prepare(size_t size) {
        if (n == size) return;
        n = size;
        a.assign(size * size, {0.0, 0.0});
        b.assign(size * size, {0.0, 0.0});
        for (auto& v : g) v.assign(size * size, 0.0f);
        fft = std::make_unique<Fft2>(size);
    }
};

Workspace& workspace(size_t n) {
    thread_local Workspace ws;
    ws.prepare(n);
    return ws;
}

// Bilinear sample; the voxel records guarantee 0 <= fx, fz < n - 1.
inline float sampleBilinear(const float* g, size_t n, float fx, float fz) {
    const size_t ix = static_cast<size_t>(fx), iz = static_cast<size_t>(fz);
    const float tx = fx - static_cast<float>(ix), tz = fz - static_cast<float>(iz);
    const float* p = g + iz * n + ix;
    return (1.0f - tz) * ((1.0f - tx) * p[0] + tx * p[1]) + tz * ((1.0f - tx) * p[n] + tx * p[n + 1]);
}

} // namespace

FluenceDoseModel::FluenceDoseModel(size_t numVoxels, size_t imageSize, double pixelMm, double bixelBoxMm)
    : m_numVoxels(numVoxels), m_n(imageSize), m_pixel(pixelMm), m_boxMm(bixelBoxMm) {
    if (imageSize < 4 || (imageSize & (imageSize - 1)) != 0)
        throw std::invalid_argument("FluenceDoseModel: image size must be a power of two");
}

void FluenceDoseModel::addKernel(size_t kernelIdx, const std::array<std::vector<double>, 3>& table,
                                 size_t tableSize) {
    m_cuda.reset();
    if (tableSize % 2 == 0 || tableSize > m_n)
        throw std::invalid_argument("FluenceDoseModel: kernel table must be odd and fit the image");
    Fft2 fft(m_n);
    Spectrum spec;
    std::vector<C> buf(m_n * m_n);
    const long half = static_cast<long>(tableSize / 2);
    for (int k = 0; k < 3; ++k) {
        std::fill(buf.begin(), buf.end(), C{0.0, 0.0});
        if (table[static_cast<size_t>(k)].size() == tableSize * tableSize) {
            for (long i = -half; i <= half; ++i) {
                for (long j = -half; j <= half; ++j) {
                    const size_t r = static_cast<size_t>((i + static_cast<long>(m_n)) % static_cast<long>(m_n));
                    const size_t c = static_cast<size_t>((j + static_cast<long>(m_n)) % static_cast<long>(m_n));
                    buf[r * m_n + c].r = table[static_cast<size_t>(k)][static_cast<size_t>(i + half) * tableSize +
                                                                          static_cast<size_t>(j + half)];
                }
            }
        }
        fft.forward(buf.data());
        auto& re = spec.re[static_cast<size_t>(k)];
        auto& im = spec.im[static_cast<size_t>(k)];
        re.resize(buf.size());
        im.resize(buf.size());
        for (size_t p = 0; p < buf.size(); ++p) { re[p] = buf[p].r; im[p] = buf[p].i; }
    }
    m_kernelSlot[kernelIdx] = m_spectra.size();
    m_spectra.push_back(std::move(spec));
}

void FluenceDoseModel::addBeam(size_t kernelIdx, std::vector<VoxelRecord> records) {
    m_cuda.reset();
    auto it = m_kernelSlot.find(kernelIdx);
    if (it == m_kernelSlot.end()) throw std::invalid_argument("FluenceDoseModel: unknown kernel");
    m_beams.push_back({it->second, std::move(records)});
}

size_t FluenceDoseModel::memoryBytes() const {
    size_t bytes = m_spectra.size() * 6 * m_n * m_n * sizeof(double);
    for (const auto& b : m_beams) bytes += b.records.capacity() * sizeof(VoxelRecord);
    return bytes;
}

bool FluenceDoseModel::enableCuda(size_t memoryBudgetBytes, bool required) {
    m_cuda.reset();
#ifdef OPTIRAD_HAS_CUDA
    try {
        m_cuda = std::make_shared<CudaFluenceBackend>(*this, memoryBudgetBytes);
        Logger::info("FluenceDoseModel: CUDA enabled, " + std::to_string(cudaMemoryBytes() / (1024 * 1024)) +
            " MiB, geometry " + (cudaGeometryResident() ? "resident" : "streamed per beam"));
        return true;
    } catch (const std::exception& error) {
        if (required) throw;
        Logger::warn("FluenceDoseModel: CUDA unavailable, using CPU: " + std::string(error.what()));
    }
#else
    (void)memoryBudgetBytes;
    if (required) throw std::runtime_error("CUDA support was not built; configure OPTIRAD_ENABLE_CUDA=ON");
#endif
    return false;
}

size_t FluenceDoseModel::cudaMemoryBytes() const {
#ifdef OPTIRAD_HAS_CUDA
    return m_cuda ? m_cuda->memoryBytes() : 0;
#else
    return 0;
#endif
}

bool FluenceDoseModel::cudaGeometryResident() const {
#ifdef OPTIRAD_HAS_CUDA
    return m_cuda && m_cuda->geometryResident();
#else
    return false;
#endif
}

void FluenceDoseModel::rasterizeBixels(const Beam& beam, const double* weights, std::vector<float>& image) const {
    image.assign(m_n * m_n, 0.0f);
    const double h = m_pixel, halfBox = 0.5 * m_boxMm;
    const long n = static_cast<long>(m_n);
    size_t bixel = 0;
    for (size_t ri = 0; ri < beam.getNumOfRays(); ++ri) {
        const auto* ray = beam.getRay(ri);
        if (!ray) continue;
        const double w = weights[bixel];
        bixel += ray->getNumOfBixels();
        if (w == 0.0) continue;

        const auto pos = ray->getRayPosBev();
        const double x0 = pos[0] - halfBox, x1 = pos[0] + halfBox;
        const double z0 = pos[2] - halfBox, z1 = pos[2] + halfBox;
        const long c0 = std::max(0L, static_cast<long>(std::floor(toPixel(x0) - 0.5)));
        const long c1 = std::min(n - 1, static_cast<long>(std::ceil(toPixel(x1) + 0.5)));
        const long r0 = std::max(0L, static_cast<long>(std::floor(toPixel(z0) - 0.5)));
        const long r1 = std::min(n - 1, static_cast<long>(std::ceil(toPixel(z1) + 0.5)));

        for (long r = r0; r <= r1; ++r) {
            const double zc = static_cast<double>(r - n / 2) * h;
            const double oz = std::max(0.0, std::min(z1, zc + 0.5 * h) - std::max(z0, zc - 0.5 * h)) / h;
            if (oz <= 0.0) continue;
            for (long c = c0; c <= c1; ++c) {
                const double xc = static_cast<double>(c - n / 2) * h;
                const double ox = std::max(0.0, std::min(x1, xc + 0.5 * h) - std::max(x0, xc - 0.5 * h)) / h;
                if (ox > 0.0) image[static_cast<size_t>(r * n + c)] += static_cast<float>(w * ox * oz);
            }
        }
    }
}

void FluenceDoseModel::gatherBixels(const Beam& beam, const std::vector<float>& image, double* out) const {
    const double h = m_pixel, halfBox = 0.5 * m_boxMm;
    const long n = static_cast<long>(m_n);
    size_t bixel = 0;
    for (size_t ri = 0; ri < beam.getNumOfRays(); ++ri) {
        const auto* ray = beam.getRay(ri);
        if (!ray) continue;
        const size_t idx = bixel;
        bixel += ray->getNumOfBixels();

        const auto pos = ray->getRayPosBev();
        const double x0 = pos[0] - halfBox, x1 = pos[0] + halfBox;
        const double z0 = pos[2] - halfBox, z1 = pos[2] + halfBox;
        const long c0 = std::max(0L, static_cast<long>(std::floor(toPixel(x0) - 0.5)));
        const long c1 = std::min(n - 1, static_cast<long>(std::ceil(toPixel(x1) + 0.5)));
        const long r0 = std::max(0L, static_cast<long>(std::floor(toPixel(z0) - 0.5)));
        const long r1 = std::min(n - 1, static_cast<long>(std::ceil(toPixel(z1) + 0.5)));

        double sum = 0.0;
        for (long r = r0; r <= r1; ++r) {
            const double zc = static_cast<double>(r - n / 2) * h;
            const double oz = std::max(0.0, std::min(z1, zc + 0.5 * h) - std::max(z0, zc - 0.5 * h)) / h;
            if (oz <= 0.0) continue;
            for (long c = c0; c <= c1; ++c) {
                const double xc = static_cast<double>(c - n / 2) * h;
                const double ox = std::max(0.0, std::min(x1, xc + 0.5 * h) - std::max(x0, xc - 0.5 * h)) / h;
                if (ox > 0.0) sum += ox * oz * static_cast<double>(image[static_cast<size_t>(r * n + c)]);
            }
        }
        out[idx] += sum;
    }
}

void FluenceDoseModel::addBeamDose(size_t beam, const std::vector<float>& fluence, double* dose) const {
    const auto& br = m_beams.at(beam);
    if (fluence.size() != m_n * m_n) throw std::invalid_argument("FluenceDoseModel: bad fluence image size");
    const Spectrum& spec = m_spectra[br.kernelSlot];
    Workspace& ws = workspace(m_n);
    const size_t total = m_n * m_n;

    for (size_t p = 0; p < total; ++p) ws.a[p] = {static_cast<double>(fluence[p]), 0.0};
    ws.fft->forward(ws.a.data());

    for (size_t k = 0; k < 3; ++k) {
        const auto& kr = spec.re[k];
        const auto& ki = spec.im[k];
        for (size_t p = 0; p < total; ++p)
            ws.b[p] = {ws.a[p].r * kr[p] - ws.a[p].i * ki[p], ws.a[p].r * ki[p] + ws.a[p].i * kr[p]};
        ws.fft->inverse(ws.b.data());
        for (size_t p = 0; p < total; ++p) ws.g[k][p] = static_cast<float>(ws.b[p].r);
    }

    for (const auto& rec : br.records) {
        const double v = rec.d[0] * sampleBilinear(ws.g[0].data(), m_n, rec.fx, rec.fz) +
                         rec.d[1] * sampleBilinear(ws.g[1].data(), m_n, rec.fx, rec.fz) +
                         rec.d[2] * sampleBilinear(ws.g[2].data(), m_n, rec.fx, rec.fz);
        dose[rec.voxel] += v;
    }
}

void FluenceDoseModel::beamGradient(size_t beam, const double* voxelGrad, std::vector<float>& fluenceGrad) const {
    const auto& br = m_beams.at(beam);
    const Spectrum& spec = m_spectra[br.kernelSlot];
    Workspace& ws = workspace(m_n);
    const size_t total = m_n * m_n;

    std::fill(ws.b.begin(), ws.b.end(), C{0.0, 0.0});
    for (size_t k = 0; k < 3; ++k) {
        std::fill(ws.a.begin(), ws.a.end(), C{0.0, 0.0});
        for (const auto& rec : br.records) {
            const double gv = voxelGrad[rec.voxel] * rec.d[k];
            if (gv == 0.0) continue;
            const size_t ix = static_cast<size_t>(rec.fx), iz = static_cast<size_t>(rec.fz);
            const double tx = rec.fx - static_cast<float>(ix), tz = rec.fz - static_cast<float>(iz);
            C* p = ws.a.data() + iz * m_n + ix;
            p[0].r += gv * (1.0 - tz) * (1.0 - tx);
            p[1].r += gv * (1.0 - tz) * tx;
            p[m_n].r += gv * tz * (1.0 - tx);
            p[m_n + 1].r += gv * tz * tx;
        }
        ws.fft->forward(ws.a.data());
        const auto& kr = spec.re[k];
        const auto& ki = spec.im[k];
        for (size_t p = 0; p < total; ++p) {
            ws.b[p].r += ws.a[p].r * kr[p] + ws.a[p].i * ki[p];
            ws.b[p].i += ws.a[p].i * kr[p] - ws.a[p].r * ki[p];
        }
    }
    ws.fft->inverse(ws.b.data());
    fluenceGrad.resize(total);
    for (size_t p = 0; p < total; ++p) fluenceGrad[p] = static_cast<float>(ws.b[p].r);
}

std::vector<double> FluenceDoseModel::computeDose(const std::vector<std::vector<float>>& fluence) const {
    if (fluence.size() != m_beams.size()) throw std::invalid_argument("FluenceDoseModel: bad beam count");
    for (const auto& image : fluence)
        if (!image.empty() && image.size() != imagePixels())
            throw std::invalid_argument("FluenceDoseModel: bad fluence image size");
#ifdef OPTIRAD_HAS_CUDA
    if (m_cuda) return m_cuda->computeDose(fluence);
#endif
    std::vector<double> dose(m_numVoxels, 0.0);
#ifdef _OPENMP
    const int nt = omp_get_max_threads();
#else
    const int nt = 1;
#endif
    std::vector<std::vector<double>> acc(static_cast<size_t>(nt));
    const long nb = static_cast<long>(std::min(m_beams.size(), fluence.size()));

#pragma omp parallel
    {
        int tid = 0;
#ifdef _OPENMP
        tid = omp_get_thread_num();
#endif
        auto& local = acc[static_cast<size_t>(tid)];
#pragma omp for schedule(dynamic, 1)
        for (long b = 0; b < nb; ++b) {
            if (fluence[static_cast<size_t>(b)].empty()) continue;
            if (local.empty()) local.assign(m_numVoxels, 0.0);
            addBeamDose(static_cast<size_t>(b), fluence[static_cast<size_t>(b)], local.data());
        }
    }
    for (const auto& local : acc)
        if (!local.empty())
            for (size_t v = 0; v < m_numVoxels; ++v) dose[v] += local[v];
    return dose;
}

std::vector<std::vector<float>> FluenceDoseModel::computeGradient(const std::vector<double>& voxelGrad) const {
    if (voxelGrad.size() != m_numVoxels) throw std::invalid_argument("FluenceDoseModel: bad voxel gradient size");
#ifdef OPTIRAD_HAS_CUDA
    if (m_cuda) return m_cuda->computeGradient(voxelGrad);
#endif
    std::vector<std::vector<float>> grads(m_beams.size());
#pragma omp parallel for schedule(dynamic, 1)
    for (long b = 0; b < static_cast<long>(m_beams.size()); ++b)
        beamGradient(static_cast<size_t>(b), voxelGrad.data(), grads[static_cast<size_t>(b)]);
    return grads;
}

} // namespace optirad
