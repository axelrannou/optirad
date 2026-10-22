#include "CudaDirectDose.hpp"
#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace optirad {
namespace {

void check(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(std::string("CUDA direct dose: ") + cudaGetErrorString(status));
}

template<class Value> struct Buffer {
    Value* data = nullptr;
    ~Buffer() { if (data) cudaFree(data); }
    void allocate(size_t count) { if (count) check(cudaMalloc(reinterpret_cast<void**>(&data), count * sizeof(Value))); }
};

struct Voxel { size_t index; double lateral, longitudinal, depthAxis, depth, distance; };
struct RayWeight { double lateral, longitudinal, weight; };
struct Kernel { double minX, minZ, dx, dz; size_t nx, nz, offset; };
struct Physics { double beta[3], mu, sad, cutoff; bool divergent; };

__device__ double interpolate(const double* data, Kernel kernel, double lateral, double longitudinal) {
    const double horizontal = (lateral - kernel.minX) / kernel.dx;
    const double vertical = (longitudinal - kernel.minZ) / kernel.dz;
    if (horizontal < 0.0 || vertical < 0.0 || horizontal > kernel.nx - 1 || vertical > kernel.nz - 1) return 0.0;
    const size_t column = min(static_cast<size_t>(horizontal), kernel.nx - 2);
    const size_t row = min(static_cast<size_t>(vertical), kernel.nz - 2);
    const double fractionX = horizontal - column, fractionZ = vertical - row;
    const double* pixel = data + kernel.offset + column * kernel.nz + row;
    return (1.0 - fractionX) * (1.0 - fractionZ) * pixel[0] + fractionX * (1.0 - fractionZ) * pixel[kernel.nz] +
        (1.0 - fractionX) * fractionZ * pixel[1] + fractionX * fractionZ * pixel[kernel.nz + 1];
}

__global__ void deposit(const Voxel* voxels, size_t count, const RayWeight* rays, size_t numRays,
                        const double* tables, const Kernel* kernels, Physics physics, double* dose) {
    const size_t index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const auto voxel = voxels[index];
    const double distanceAxis = physics.sad + voxel.depthAxis;
    if (distanceAxis < 1.0 || voxel.depth <= 0.0 || voxel.distance < 1.0) return;
    double depthDose[3];
    for (int component = 0; component < 3; ++component) {
        const double beta = physics.beta[component], denominator = beta - physics.mu;
        depthDose[component] = fabs(denominator) > 1e-12 ? beta / denominator *
            (exp(-physics.mu * voxel.depth) - exp(-beta * voxel.depth)) :
            physics.mu * voxel.depth * exp(-physics.mu * voxel.depth);
    }
    const double inverseSquare = physics.sad * physics.sad / (voxel.distance * voxel.distance);
    double sum = 0.0;
    for (size_t ray = 0; ray < numRays; ++ray) {
        if (rays[ray].weight <= 0.0) continue;
        const double lateral = physics.divergent ? voxel.lateral * physics.sad / distanceAxis - rays[ray].lateral :
            (voxel.lateral - rays[ray].lateral) * physics.sad / distanceAxis;
        const double longitudinal = physics.divergent ? voxel.longitudinal * physics.sad / distanceAxis - rays[ray].longitudinal :
            (voxel.longitudinal - rays[ray].longitudinal) * physics.sad / distanceAxis;
        if (sqrt(lateral * lateral + longitudinal * longitudinal) > physics.cutoff) continue;
        double value = 0.0;
        for (int component = 0; component < 3; ++component)
            value += depthDose[component] * interpolate(tables, kernels[component], lateral, longitudinal);
        sum += rays[ray].weight * fmax(0.0, value * inverseSquare);
    }
    dose[voxel.index] += sum;
}

} // namespace

struct CudaDirectDose::Impl {
    size_t count = 0;
    int device = 0;
    Buffer<double> dose;
    ~Impl() { cudaSetDevice(device); }
};

CudaDirectDose::CudaDirectDose(size_t numVoxels) : m_impl(std::make_unique<Impl>()) {
    check(cudaGetDevice(&m_impl->device));
    m_impl->count = numVoxels;
    m_impl->dose.allocate(numVoxels);
    check(cudaMemset(m_impl->dose.data, 0, numVoxels * sizeof(double)));
}
CudaDirectDose::~CudaDirectDose() = default;

void CudaDirectDose::accumulate(const Beam& beam, const std::vector<size_t>& indices,
    const std::vector<Vec3>& coordinates, const std::vector<double>& depths, const std::vector<double>& distances,
    const GridInterpolant2D* kernels, const std::vector<double>& weights, const std::array<double, 3>& betas,
    double mu, double sad, double cutoff, bool divergent) {
    check(cudaSetDevice(m_impl->device));
    if (indices.empty()) return;
    std::vector<Voxel> voxels;
    voxels.reserve(indices.size());
    for (size_t index = 0; index < indices.size(); ++index)
        voxels.push_back({indices[index], coordinates[index][0], coordinates[index][2], coordinates[index][1], depths[index], distances[index]});
    std::vector<RayWeight> rays;
    for (size_t index = 0; index < beam.getNumOfRays(); ++index) {
        const auto position = beam.getRay(index)->getRayPosBev();
        rays.push_back({position[0], position[2], weights[index]});
    }
    if (rays.empty()) return;
    std::vector<double> tables;
    Kernel descriptions[3];
    for (size_t component = 0; component < 3; ++component) {
        const auto& kernel = kernels[component];
        descriptions[component] = {kernel.m_xMin, kernel.m_yMin, kernel.m_dx, kernel.m_dy, kernel.m_nx, kernel.m_ny, tables.size()};
        tables.insert(tables.end(), kernel.m_data.begin(), kernel.m_data.end());
    }
    Buffer<Voxel> deviceVoxels;
    Buffer<RayWeight> deviceRays;
    Buffer<double> deviceTables;
    Buffer<Kernel> deviceKernels;
    deviceVoxels.allocate(voxels.size()); deviceRays.allocate(rays.size());
    deviceTables.allocate(tables.size()); deviceKernels.allocate(3);
    check(cudaMemcpy(deviceVoxels.data, voxels.data(), voxels.size() * sizeof(Voxel), cudaMemcpyHostToDevice));
    check(cudaMemcpy(deviceRays.data, rays.data(), rays.size() * sizeof(RayWeight), cudaMemcpyHostToDevice));
    check(cudaMemcpy(deviceTables.data, tables.data(), tables.size() * sizeof(double), cudaMemcpyHostToDevice));
    check(cudaMemcpy(deviceKernels.data, descriptions, sizeof(descriptions), cudaMemcpyHostToDevice));
    Physics physics{{betas[0], betas[1], betas[2]}, mu, sad, cutoff, divergent};
    deposit<<<static_cast<unsigned>((voxels.size() + 127) / 128), 128>>>(deviceVoxels.data, voxels.size(),
        deviceRays.data, rays.size(), deviceTables.data, deviceKernels.data, physics, m_impl->dose.data);
    check(cudaGetLastError());
    check(cudaDeviceSynchronize());
}

std::vector<double> CudaDirectDose::result() {
    check(cudaSetDevice(m_impl->device));
    std::vector<double> output(m_impl->count);
    check(cudaMemcpy(output.data(), m_impl->dose.data, output.size() * sizeof(double), cudaMemcpyDeviceToHost));
    return output;
}

} // namespace optirad