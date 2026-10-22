#include "CudaFluenceBackend.hpp"
#include "FluenceDoseModel.hpp"
#include "core/Stf.hpp"
#include "IDoseOperator.hpp"
#include <thrust/device_ptr.h>
#include <thrust/sort.h>

#include <cuda_runtime.h>
#include <cufft.h>
#include <algorithm>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>

namespace optirad {
namespace {

void checkCuda(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(std::string("CUDA: ") + cudaGetErrorString(status));
}

void checkFft(cufftResult status) {
    if (status != CUFFT_SUCCESS) throw std::runtime_error("cuFFT error " + std::to_string(status));
}

template<class Value> struct DeviceBuffer {
    Value* data = nullptr;
    DeviceBuffer() = default;
    DeviceBuffer(const DeviceBuffer&) = delete;
    DeviceBuffer& operator=(const DeviceBuffer&) = delete;
    ~DeviceBuffer() { if (data) cudaFree(data); }
    void allocate(size_t count) {
        if (count) checkCuda(cudaMalloc(reinterpret_cast<void**>(&data), count * sizeof(Value)));
    }
};

struct FftPlan {
    cufftHandle handle = 0;
    ~FftPlan() { if (handle) cufftDestroy(handle); }
};

struct BixelBox { double left, right, bottom, top; };
struct PixelLink { size_t bixel; double overlapX, overlapZ; };

__global__ void rasterizeBixels(const size_t* offsets, const PixelLink* links, const double* weights,
                               size_t pixels, float* image) {
    const size_t pixel = blockIdx.x * blockDim.x + threadIdx.x;
    if (pixel >= pixels) return;
    float value = 0.0f;
    for (size_t index = offsets[pixel]; index < offsets[pixel + 1]; ++index) {
        const auto link = links[index];
        value += static_cast<float>(weights[link.bixel] * link.overlapX * link.overlapZ);
    }
    image[pixel] = value;
}

__global__ void gatherBixels(const BixelBox* boxes, size_t count, size_t width, double pixelMm,
                            const float* image, double* gradient) {
    const size_t bixel = blockIdx.x * blockDim.x + threadIdx.x;
    if (bixel >= count) return;
    const auto box = boxes[bixel];
    const int half = static_cast<int>(width / 2);
    const int columnStart = max(0, static_cast<int>(floor(box.left / pixelMm + half - 0.5)));
    const int columnStop = min(static_cast<int>(width) - 1, static_cast<int>(ceil(box.right / pixelMm + half + 0.5)));
    const int rowStart = max(0, static_cast<int>(floor(box.bottom / pixelMm + half - 0.5)));
    const int rowStop = min(static_cast<int>(width) - 1, static_cast<int>(ceil(box.top / pixelMm + half + 0.5)));
    double value = 0.0;
    for (int row = rowStart; row <= rowStop; ++row) {
        const double centerZ = (row - half) * pixelMm;
        const double overlapZ = fmax(0.0, fmin(box.top, centerZ + pixelMm * 0.5) - fmax(box.bottom, centerZ - pixelMm * 0.5)) / pixelMm;
        for (int column = columnStart; column <= columnStop; ++column) {
            const double centerX = (column - half) * pixelMm;
            const double overlapX = fmax(0.0, fmin(box.right, centerX + pixelMm * 0.5) - fmax(box.left, centerX - pixelMm * 0.5)) / pixelMm;
            value += overlapX * overlapZ * image[row * width + column];
        }
    }
    gradient[bixel] = value;
}

__global__ void loadImage(const float* image, cufftDoubleComplex* spectrum, size_t pixels) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < pixels) spectrum[index] = {static_cast<double>(image[index]), 0.0};
}

__global__ void multiplySpectrum(const cufftDoubleComplex* input, const cufftDoubleComplex* kernel,
                                cufftDoubleComplex* output, size_t pixels, bool adjoint, bool accumulate) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= pixels) return;
    const auto value = input[index];
    auto factor = kernel[index];
    if (adjoint) factor.y = -factor.y;
    const cufftDoubleComplex product{value.x * factor.x - value.y * factor.y,
                                     value.x * factor.y + value.y * factor.x};
    if (accumulate) {
        output[index].x += product.x;
        output[index].y += product.y;
    } else output[index] = product;
}

__global__ void storeImage(const cufftDoubleComplex* spectrum, float* image, size_t pixels) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < pixels) image[index] = static_cast<float>(spectrum[index].x / static_cast<double>(pixels));
}

__device__ float sample(const float* image, size_t width, float lateral, float longitudinal) {
    const size_t column = static_cast<size_t>(lateral), row = static_cast<size_t>(longitudinal);
    const float horizontal = lateral - static_cast<float>(column);
    const float vertical = longitudinal - static_cast<float>(row);
    const float* pixel = image + row * width + column;
    return (1.0f - vertical) * ((1.0f - horizontal) * pixel[0] + horizontal * pixel[1]) +
           vertical * ((1.0f - horizontal) * pixel[width] + horizontal * pixel[width + 1]);
}

__global__ void gatherDose(const FluenceDoseModel::VoxelRecord* records, size_t count,
                          const float* images, size_t width, double* dose) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const auto record = records[index];
    const size_t pixels = width * width;
    const double value = record.d[0] * sample(images, width, record.fx, record.fz) +
        record.d[1] * sample(images + pixels, width, record.fx, record.fz) +
        record.d[2] * sample(images + 2 * pixels, width, record.fx, record.fz);
    atomicAdd(dose + record.voxel, value);
}

__global__ void scatterGradient(const FluenceDoseModel::VoxelRecord* records, size_t count,
                               const double* gradient, size_t width, int component, cufftDoubleComplex* image) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const auto record = records[index];
    const double value = gradient[record.voxel] * record.d[component];
    const size_t column = static_cast<size_t>(record.fx), row = static_cast<size_t>(record.fz);
    const double horizontal = record.fx - static_cast<float>(column);
    const double vertical = record.fz - static_cast<float>(row);
    auto* pixel = image + row * width + column;
    atomicAdd(&pixel[0].x, value * (1.0 - vertical) * (1.0 - horizontal));
    atomicAdd(&pixel[1].x, value * (1.0 - vertical) * horizontal);
    atomicAdd(&pixel[width].x, value * vertical * (1.0 - horizontal));
    atomicAdd(&pixel[width + 1].x, value * vertical * horizontal);
}

unsigned blocks(size_t count) { return static_cast<unsigned>((count + 255) / 256); }

struct DoseRank { double dose; size_t voxel; };
struct DoseRankLess {
    __host__ __device__ bool operator()(const DoseRank& first, const DoseRank& second) const {
        return first.dose < second.dose || (first.dose == second.dose && first.voxel < second.voxel);
    }
};

__global__ void rankDoses(const size_t* indices, size_t count, const double* dose, DoseRank* ranked) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index < count) ranked[index] = {dose[indices[index]], indices[index]};
}

__global__ void doseObjective(const size_t* indices, const DoseRank* ranked, size_t count, size_t begin,
                             size_t end, int type, double threshold, double scale, const double* dose,
                             double* gradient, double* value) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count || index < begin || index >= end) return;
    const size_t voxel = type >= 3 ? ranked[index].voxel : indices[index];
    const double difference = dose[voxel] - threshold;
    if ((type == 1 || type == 4) && difference <= 0.0) return;
    if ((type == 2 || type == 3) && difference >= 0.0) return;
    atomicAdd(value, scale * difference * difference);
    atomicAdd(gradient + voxel, 2.0 * scale * difference);
}

__global__ void hotspotObjective(const double* dose, size_t count, double prescription, double threshold,
                                double penalty, double* gradient, double* value, bool consistent) {
    const size_t voxel = blockIdx.x * blockDim.x + threadIdx.x;
    if (voxel >= count) return;
    const double excess = dose[voxel] - prescription * threshold;
    if (excess <= 0.0) return;
    const double normalized = excess / prescription;
    atomicAdd(value, penalty * normalized * normalized * (1.0 + 10.0 * normalized) * excess * excess);
    gradient[voxel] += consistent ? penalty * excess * normalized * normalized * (4.0 + 50.0 * normalized) :
        2.0 * penalty * normalized * (1.0 + 15.0 * normalized) * excess;
}

} // namespace

struct CudaFluenceBackend::Impl {
    size_t width = 0, pixels = 0, voxels = 0, bytes = 0;
    int device = 0;
    bool resident = false;
    size_t budget = 0;
    FftPlan plan;
    DeviceBuffer<unsigned char> fftWorkspace;
    DeviceBuffer<cufftDoubleComplex> kernels, input, work;
    DeviceBuffer<float> images;
    DeviceBuffer<double> dose;
    DeviceBuffer<FluenceDoseModel::VoxelRecord> geometry;
    struct BeamInfo {
        size_t slot = 0, offset = 0, count = 0;
        std::vector<FluenceDoseModel::VoxelRecord> hostRecords;
    };
    std::vector<BeamInfo> beams;
    std::mutex mutex;

    ~Impl() { cudaSetDevice(device); }

    const FluenceDoseModel::VoxelRecord* records(const BeamInfo& beam) {
        if (!resident && beam.count)
            checkCuda(cudaMemcpy(geometry.data, beam.hostRecords.data(),
                beam.count * sizeof(FluenceDoseModel::VoxelRecord), cudaMemcpyHostToDevice));
        return geometry.data ? geometry.data + (resident ? beam.offset : 0) : nullptr;
    }

    void fft(cufftDoubleComplex* buffer, int direction) {
        checkFft(cufftExecZ2Z(plan.handle, buffer, buffer, direction));
    }
    void bixelForward(BixelLayout& layout, const std::vector<double>& weights);
    void bixelAdjoint(BixelLayout& layout);

    void forward(const BeamInfo& beam) {
        loadImage<<<blocks(pixels), 256>>>(images.data, input.data, pixels);
        fft(input.data, CUFFT_FORWARD);
        for (size_t component = 0; component < 3; ++component) {
            multiplySpectrum<<<blocks(pixels), 256>>>(input.data,
                kernels.data + (beam.slot * 3 + component) * pixels, work.data, pixels, false, false);
            fft(work.data, CUFFT_INVERSE);
            storeImage<<<blocks(pixels), 256>>>(work.data, images.data + component * pixels, pixels);
        }
        gatherDose<<<blocks(beam.count), 256>>>(records(beam), beam.count, images.data, width, dose.data);
        checkCuda(cudaGetLastError());
    }

    void adjoint(const BeamInfo& beam) {
        const auto* recordData = records(beam);
        checkCuda(cudaMemset(input.data, 0, pixels * sizeof(cufftDoubleComplex)));
        for (int component = 0; component < 3; ++component) {
            checkCuda(cudaMemset(work.data, 0, pixels * sizeof(cufftDoubleComplex)));
            scatterGradient<<<blocks(beam.count), 256>>>(recordData, beam.count, dose.data,
                width, component, work.data);
            fft(work.data, CUFFT_FORWARD);
            multiplySpectrum<<<blocks(pixels), 256>>>(work.data,
                kernels.data + (beam.slot * 3 + component) * pixels, input.data, pixels, true, true);
        }
        fft(input.data, CUFFT_INVERSE);
        storeImage<<<blocks(pixels), 256>>>(input.data, images.data, pixels);
        checkCuda(cudaGetLastError());
    }
};

struct CudaFluenceBackend::BixelLayout {
    int device = 0;
    size_t count = 0;
    double pixelMm = 0.0;
    std::vector<size_t> offsets;
    DeviceBuffer<BixelBox> boxes;
    DeviceBuffer<size_t> pixelOffsets;
    DeviceBuffer<PixelLink> pixelLinks;
    DeviceBuffer<double> values;
    size_t remainingBudget = 0;
    struct Objectives {
        std::vector<std::vector<size_t>> indices;
        std::vector<size_t> offsets;
        DeviceBuffer<size_t> deviceIndices;
        DeviceBuffer<DoseRank> ranked;
        DeviceBuffer<double> gradient, value;
    };
    std::unique_ptr<Objectives> objectives;
    ~BixelLayout() { cudaSetDevice(device); }
};

CudaFluenceBackend::CudaFluenceBackend(const FluenceDoseModel& model, size_t memoryBudgetBytes)
    : m_impl(std::make_unique<Impl>()) {
    auto& state = *m_impl;
    checkCuda(cudaGetDevice(&state.device));
    state.width = model.m_n;
    state.pixels = model.imagePixels();
    state.voxels = model.numVoxels();
    if (!state.voxels || state.width > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("CUDA: unsupported dose or image dimensions");
    checkFft(cufftCreate(&state.plan.handle));
    checkFft(cufftSetAutoAllocation(state.plan.handle, 0));
    size_t workspaceBytes = 0;
    checkFft(cufftMakePlan2d(state.plan.handle, static_cast<int>(state.width),
        static_cast<int>(state.width), CUFFT_Z2Z, &workspaceBytes));
    size_t freeBytes = 0, totalBytes = 0;
    checkCuda(cudaMemGetInfo(&freeBytes, &totalBytes));
    const size_t budget = memoryBudgetBytes ? std::min(memoryBudgetBytes, freeBytes * 8 / 10) : freeBytes * 7 / 10;
    state.budget = budget;
    const size_t kernelElements = model.m_spectra.size() * 3 * state.pixels;
    const size_t fixedBytes = workspaceBytes + (kernelElements + 2 * state.pixels) * sizeof(cufftDoubleComplex) +
        3 * state.pixels * sizeof(float) + state.voxels * sizeof(double);
    size_t totalRecords = 0, maxRecords = 0;
    for (const auto& beam : model.m_beams) {
        totalRecords += beam.records.size();
        maxRecords = std::max(maxRecords, beam.records.size());
    }
    if (fixedBytes > budget || maxRecords > (budget - fixedBytes) / sizeof(FluenceDoseModel::VoxelRecord))
        throw std::runtime_error("CUDA memory budget cannot hold workspace and one beam; increase budget or use CPU");
    state.resident = totalRecords <= (budget - fixedBytes) / sizeof(FluenceDoseModel::VoxelRecord);
    state.bytes = fixedBytes + (state.resident ? totalRecords : maxRecords) * sizeof(FluenceDoseModel::VoxelRecord);
    state.fftWorkspace.allocate(workspaceBytes);
    checkFft(cufftSetWorkArea(state.plan.handle, state.fftWorkspace.data));
    state.kernels.allocate(kernelElements);
    state.input.allocate(state.pixels);
    state.work.allocate(state.pixels);
    state.images.allocate(3 * state.pixels);
    state.dose.allocate(state.voxels);
    state.geometry.allocate(state.resident ? totalRecords : maxRecords);
    std::vector<cufftDoubleComplex> spectrum(state.pixels);
    for (size_t slot = 0; slot < model.m_spectra.size(); ++slot) {
        for (size_t component = 0; component < 3; ++component) {
            for (size_t index = 0; index < state.pixels; ++index)
                spectrum[index] = {model.m_spectra[slot].re[component][index], model.m_spectra[slot].im[component][index]};
            checkCuda(cudaMemcpy(state.kernels.data + (slot * 3 + component) * state.pixels,
                spectrum.data(), state.pixels * sizeof(cufftDoubleComplex), cudaMemcpyHostToDevice));
        }
    }
    size_t offset = 0;
    for (const auto& beam : model.m_beams) {
        Impl::BeamInfo info;
        info.slot = beam.kernelSlot;
        info.offset = offset;
        info.count = beam.records.size();
        if (state.resident && info.count)
            checkCuda(cudaMemcpy(state.geometry.data + offset, beam.records.data(),
                info.count * sizeof(FluenceDoseModel::VoxelRecord), cudaMemcpyHostToDevice));
        if (!state.resident) info.hostRecords = beam.records;
        state.beams.push_back(std::move(info));
        offset += beam.records.size();
    }
}

CudaFluenceBackend::~CudaFluenceBackend() = default;
size_t CudaFluenceBackend::memoryBytes() const { return m_impl->bytes; }
bool CudaFluenceBackend::geometryResident() const { return m_impl->resident; }

std::shared_ptr<CudaFluenceBackend::BixelLayout> CudaFluenceBackend::createBixelLayout(
    const Stf& stf, double pixelMm, double boxMm) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    auto layout = std::make_shared<BixelLayout>();
    layout->device = state.device;
    layout->pixelMm = pixelMm;
    std::vector<BixelBox> boxes;
    layout->offsets.push_back(0);
    for (const auto& beam : stf.getBeams()) {
        for (const auto& ray : beam.getRays()) {
            if (ray.getNumOfBixels() != 1) throw std::invalid_argument("CUDA bixel mapping requires photon rays");
            const auto position = ray.getRayPosBev();
            boxes.push_back({position[0] - boxMm * 0.5, position[0] + boxMm * 0.5,
                position[2] - boxMm * 0.5, position[2] + boxMm * 0.5});
        }
        layout->offsets.push_back(boxes.size());
    }
    layout->count = boxes.size();
    std::vector<size_t> pixelOffsets;
    std::vector<PixelLink> pixelLinks;
    for (size_t beam = 0; beam < state.beams.size(); ++beam) {
        std::vector<std::vector<PixelLink>> links(state.pixels);
        for (size_t bixel = layout->offsets[beam]; bixel < layout->offsets[beam + 1]; ++bixel) {
            const auto box = boxes[bixel];
            const long half = static_cast<long>(state.width / 2);
            const long columnStart = std::max(0L, static_cast<long>(std::floor(box.left / pixelMm + half - 0.5)));
            const long columnStop = std::min(static_cast<long>(state.width) - 1, static_cast<long>(std::ceil(box.right / pixelMm + half + 0.5)));
            const long rowStart = std::max(0L, static_cast<long>(std::floor(box.bottom / pixelMm + half - 0.5)));
            const long rowStop = std::min(static_cast<long>(state.width) - 1, static_cast<long>(std::ceil(box.top / pixelMm + half + 0.5)));
            for (long row = rowStart; row <= rowStop; ++row) {
                const double centerZ = (row - half) * pixelMm;
                const double overlapZ = std::max(0.0, std::min(box.top, centerZ + pixelMm * 0.5) - std::max(box.bottom, centerZ - pixelMm * 0.5)) / pixelMm;
                for (long column = columnStart; column <= columnStop; ++column) {
                    const double centerX = (column - half) * pixelMm;
                    const double overlapX = std::max(0.0, std::min(box.right, centerX + pixelMm * 0.5) - std::max(box.left, centerX - pixelMm * 0.5)) / pixelMm;
                    if (overlapX > 0.0 && overlapZ > 0.0)
                        links[row * state.width + column].push_back({bixel, overlapX, overlapZ});
                }
            }
        }
        for (const auto& contributors : links) {
            pixelOffsets.push_back(pixelLinks.size());
            pixelLinks.insert(pixelLinks.end(), contributors.begin(), contributors.end());
        }
        pixelOffsets.push_back(pixelLinks.size());
    }
    const size_t layoutBytes = layout->count * (sizeof(BixelBox) + sizeof(double)) +
        pixelOffsets.size() * sizeof(size_t) + pixelLinks.size() * sizeof(PixelLink);
    if (layoutBytes > state.budget - state.bytes)
        throw std::runtime_error("CUDA memory budget cannot hold bixel layout; increase budget or use CPU");
    layout->remainingBudget = state.budget - state.bytes - layoutBytes;
    layout->boxes.allocate(layout->count);
    layout->values.allocate(layout->count);
    layout->pixelOffsets.allocate(pixelOffsets.size());
    layout->pixelLinks.allocate(pixelLinks.size());
    if (!pixelOffsets.empty()) checkCuda(cudaMemcpy(layout->pixelOffsets.data, pixelOffsets.data(),
        pixelOffsets.size() * sizeof(size_t), cudaMemcpyHostToDevice));
    if (!pixelLinks.empty()) checkCuda(cudaMemcpy(layout->pixelLinks.data, pixelLinks.data(),
        pixelLinks.size() * sizeof(PixelLink), cudaMemcpyHostToDevice));
    if (layout->count) checkCuda(cudaMemcpy(layout->boxes.data, boxes.data(),
        layout->count * sizeof(BixelBox), cudaMemcpyHostToDevice));
    return layout;
}

void CudaFluenceBackend::Impl::bixelForward(BixelLayout& layout, const std::vector<double>& weights) {
    auto& state = *this;
    if (weights.size() != layout.count) throw std::invalid_argument("CUDA: bad bixel weight count");
    if (layout.count) checkCuda(cudaMemcpy(layout.values.data, weights.data(), layout.count * sizeof(double), cudaMemcpyHostToDevice));
    checkCuda(cudaMemset(state.dose.data, 0, state.voxels * sizeof(double)));
    for (size_t beam = 0; beam < state.beams.size(); ++beam) {
        const size_t begin = layout.offsets[beam], count = layout.offsets[beam + 1] - begin;
        if (!count || !state.beams[beam].count) continue;
        if (std::none_of(weights.begin() + begin, weights.begin() + begin + count,
            [](double weight) { return weight != 0.0; })) continue;
        rasterizeBixels<<<blocks(state.pixels), 256>>>(layout.pixelOffsets.data + beam * (state.pixels + 1),
            layout.pixelLinks.data, layout.values.data, state.pixels, state.images.data);
        state.forward(state.beams[beam]);
    }
}

std::vector<double> CudaFluenceBackend::computeBixelDose(BixelLayout& layout, const std::vector<double>& weights) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    state.bixelForward(layout, weights);
    std::vector<double> result(state.voxels);
    checkCuda(cudaMemcpy(result.data(), state.dose.data, state.voxels * sizeof(double), cudaMemcpyDeviceToHost));
    return result;
}

void CudaFluenceBackend::Impl::bixelAdjoint(BixelLayout& layout) {
    auto& state = *this;
    if (layout.count) checkCuda(cudaMemset(layout.values.data, 0, layout.count * sizeof(double)));
    for (size_t beam = 0; beam < state.beams.size(); ++beam) {
        const size_t begin = layout.offsets[beam], count = layout.offsets[beam + 1] - begin;
        if (!count || !state.beams[beam].count) continue;
        state.adjoint(state.beams[beam]);
        gatherBixels<<<blocks(count), 256>>>(layout.boxes.data + begin, count, state.width,
            layout.pixelMm, state.images.data, layout.values.data + begin);
        checkCuda(cudaGetLastError());
    }
}

void CudaFluenceBackend::accumulateBixelGradient(BixelLayout& layout, const std::vector<double>& voxelGradient,
                                               std::vector<double>& gradient) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    if (voxelGradient.size() != state.voxels || gradient.size() != layout.count)
        throw std::invalid_argument("CUDA: bad bixel gradient dimensions");
    checkCuda(cudaMemcpy(state.dose.data, voxelGradient.data(), state.voxels * sizeof(double), cudaMemcpyHostToDevice));
    state.bixelAdjoint(layout);
    std::vector<double> contribution(layout.count);
    if (layout.count) checkCuda(cudaMemcpy(contribution.data(), layout.values.data,
        layout.count * sizeof(double), cudaMemcpyDeviceToHost));
    for (size_t bixel = 0; bixel < layout.count; ++bixel) gradient[bixel] += contribution[bixel];
}

bool CudaFluenceBackend::evaluateObjectives(BixelLayout& layout, const std::vector<double>& weights,
    const std::vector<DoseObjectiveSpec>& objectives, double prescription, double threshold, double penalty,
    double& value, std::vector<double>& gradient, bool consistentHotspotGradient, bool needGradient) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    bool cached = layout.objectives && layout.objectives->indices.size() == objectives.size();
    size_t totalIndices = 0, maxIndices = 0;
    for (size_t index = 0; index < objectives.size(); ++index) {
        const auto& objective = objectives[index];
        if (objective.type < 0 || objective.type > 4 || objective.volumeFraction < 0.0 || objective.volumeFraction > 1.0)
            return false;
        if (cached && layout.objectives->indices[index] != objective.indices) cached = false;
        totalIndices += objective.indices.size();
        maxIndices = std::max(maxIndices, objective.indices.size());
    }
    if (!cached) {
        const size_t bytes = totalIndices * sizeof(size_t) + maxIndices * sizeof(DoseRank) * 3 +
            (state.voxels + 1) * sizeof(double);
        if (bytes > layout.remainingBudget) return false;
        layout.objectives.reset();
        auto workspace = std::make_unique<BixelLayout::Objectives>();
        std::vector<size_t> indices;
        workspace->offsets.push_back(0);
        for (const auto& objective : objectives) {
            workspace->indices.push_back(objective.indices);
            std::vector<size_t> valid;
            for (size_t voxel : objective.indices) if (voxel < state.voxels) valid.push_back(voxel);
            auto unique = valid;
            std::sort(unique.begin(), unique.end());
            if (std::adjacent_find(unique.begin(), unique.end()) != unique.end()) return false;
            indices.insert(indices.end(), valid.begin(), valid.end());
            workspace->offsets.push_back(indices.size());
        }
        workspace->deviceIndices.allocate(indices.size());
        workspace->ranked.allocate(maxIndices);
        workspace->gradient.allocate(state.voxels);
        workspace->value.allocate(1);
        if (!indices.empty()) checkCuda(cudaMemcpy(workspace->deviceIndices.data, indices.data(),
            indices.size() * sizeof(size_t), cudaMemcpyHostToDevice));
        layout.objectives = std::move(workspace);
    }
    auto& workspace = *layout.objectives;
    state.bixelForward(layout, weights);
    checkCuda(cudaMemset(workspace.gradient.data, 0, state.voxels * sizeof(double)));
    checkCuda(cudaMemset(workspace.value.data, 0, sizeof(double)));
    for (size_t index = 0; index < objectives.size(); ++index) {
        const auto& objective = objectives[index];
        const size_t offset = workspace.offsets[index], count = workspace.offsets[index + 1] - offset;
        if (!count) continue;
        const size_t* indices = workspace.deviceIndices.data + offset;
        size_t begin = 0, end = count;
        double scale = objective.weight;
        if (objective.type >= 3) {
            rankDoses<<<blocks(count), 256>>>(indices, count, state.dose.data, workspace.ranked.data);
            auto ranked = thrust::device_pointer_cast(workspace.ranked.data);
            thrust::sort(ranked, ranked + count, DoseRankLess{});
            if (objective.type == 3) begin = static_cast<size_t>(std::floor((1.0 - objective.volumeFraction) * count));
            else end = count - std::min(count, static_cast<size_t>(std::ceil(objective.volumeFraction * count)));
        } else scale /= objective.indices.size();
        doseObjective<<<blocks(count), 256>>>(indices, workspace.ranked.data, count, begin, end,
            objective.type, objective.threshold, scale, state.dose.data, workspace.gradient.data, workspace.value.data);
    }
    if (prescription > 0.0 && penalty > 0.0)
        hotspotObjective<<<blocks(state.voxels), 256>>>(state.dose.data, state.voxels, prescription, threshold,
            penalty, workspace.gradient.data, workspace.value.data, consistentHotspotGradient);
    checkCuda(cudaGetLastError());
    checkCuda(cudaMemcpy(&value, workspace.value.data, sizeof(double), cudaMemcpyDeviceToHost));
    if (!needGradient) { gradient.clear(); return true; }
    checkCuda(cudaMemcpy(state.dose.data, workspace.gradient.data, state.voxels * sizeof(double), cudaMemcpyDeviceToDevice));
    state.bixelAdjoint(layout);
    gradient.resize(layout.count);
    if (layout.count) checkCuda(cudaMemcpy(gradient.data(), layout.values.data, layout.count * sizeof(double), cudaMemcpyDeviceToHost));
    return true;
}

std::vector<double> CudaFluenceBackend::computeDose(const std::vector<std::vector<float>>& fluence) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    checkCuda(cudaMemset(state.dose.data, 0, state.voxels * sizeof(double)));
    for (size_t beamIndex = 0; beamIndex < state.beams.size(); ++beamIndex) {
        if (fluence[beamIndex].empty()) continue;
        const auto& beam = state.beams[beamIndex];
        if (!beam.count) continue;
        const auto* records = state.records(beam);
        checkCuda(cudaMemcpy(state.images.data, fluence[beamIndex].data(), state.pixels * sizeof(float), cudaMemcpyHostToDevice));
        loadImage<<<blocks(state.pixels), 256>>>(state.images.data, state.input.data, state.pixels);
        state.fft(state.input.data, CUFFT_FORWARD);
        for (size_t component = 0; component < 3; ++component) {
            multiplySpectrum<<<blocks(state.pixels), 256>>>(state.input.data,
                state.kernels.data + (beam.slot * 3 + component) * state.pixels, state.work.data, state.pixels, false, false);
            state.fft(state.work.data, CUFFT_INVERSE);
            storeImage<<<blocks(state.pixels), 256>>>(state.work.data, state.images.data + component * state.pixels, state.pixels);
        }
        gatherDose<<<blocks(beam.count), 256>>>(records, beam.count, state.images.data, state.width, state.dose.data);
        checkCuda(cudaGetLastError());
    }
    std::vector<double> result(state.voxels);
    checkCuda(cudaMemcpy(result.data(), state.dose.data, state.voxels * sizeof(double), cudaMemcpyDeviceToHost));
    return result;
}

std::vector<std::vector<float>> CudaFluenceBackend::computeGradient(const std::vector<double>& voxelGradient) {
    auto& state = *m_impl;
    std::lock_guard<std::mutex> guard(state.mutex);
    checkCuda(cudaSetDevice(state.device));
    checkCuda(cudaMemcpy(state.dose.data, voxelGradient.data(), state.voxels * sizeof(double), cudaMemcpyHostToDevice));
    std::vector<std::vector<float>> result(state.beams.size(), std::vector<float>(state.pixels, 0.0f));
    for (size_t beamIndex = 0; beamIndex < state.beams.size(); ++beamIndex) {
        const auto& beam = state.beams[beamIndex];
        if (!beam.count) continue;
        const auto* records = state.records(beam);
        checkCuda(cudaMemset(state.input.data, 0, state.pixels * sizeof(cufftDoubleComplex)));
        for (int component = 0; component < 3; ++component) {
            checkCuda(cudaMemset(state.work.data, 0, state.pixels * sizeof(cufftDoubleComplex)));
            scatterGradient<<<blocks(beam.count), 256>>>(records, beam.count, state.dose.data,
                state.width, component, state.work.data);
            state.fft(state.work.data, CUFFT_FORWARD);
            multiplySpectrum<<<blocks(state.pixels), 256>>>(state.work.data,
                state.kernels.data + (beam.slot * 3 + component) * state.pixels, state.input.data, state.pixels, true, true);
        }
        state.fft(state.input.data, CUFFT_INVERSE);
        storeImage<<<blocks(state.pixels), 256>>>(state.input.data, state.images.data, state.pixels);
        checkCuda(cudaGetLastError());
        checkCuda(cudaMemcpy(result[beamIndex].data(), state.images.data, state.pixels * sizeof(float), cudaMemcpyDeviceToHost));
    }
    return result;
}

} // namespace optirad