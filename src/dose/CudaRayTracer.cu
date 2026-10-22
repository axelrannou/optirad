#include "CudaRayTracer.hpp"
#include <cuda_runtime.h>
#include <math_constants.h>
#include <stdexcept>
#include <string>

namespace optirad {
namespace {

void check(cudaError_t status) {
    if (status != cudaSuccess) throw std::runtime_error(std::string("CUDA ray tracing: ") + cudaGetErrorString(status));
}

template<class Value> struct Buffer {
    Value* data = nullptr;
    ~Buffer() { if (data) cudaFree(data); }
    void allocate(size_t count) { if (count) check(cudaMalloc(reinterpret_cast<void**>(&data), count * sizeof(Value))); }
};

struct TraceGrid {
    int dimensions[3];
    double spacing[3], origin[3], inverse[9];
};

struct Point { double coordinate[3]; };

__global__ void traceDepths(TraceGrid grid, Point source, const Point* targets, size_t count,
                            const double* density, double* output) {
    const size_t ray = blockIdx.x * blockDim.x + threadIdx.x;
    if (ray >= count) return;
    double direction[3], distanceSquared = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        direction[axis] = targets[ray].coordinate[axis] - source.coordinate[axis];
        distanceSquared += direction[axis] * direction[axis];
    }
    output[ray] = 0.0;
    const double distance = sqrt(distanceSquared);
    if (distance < 1e-10) return;
    double alphaMin = 0.0, alphaMax = 1.0;
    int currentPlane[3] = {0, 0, 0}, lastPlane[3] = {-1, -1, -1}, increment[3] = {1, 1, 1};
    bool hasPlanes = false;
    for (int axis = 0; axis < 3; ++axis) {
        const double delta = direction[axis];
        if (fabs(delta) < 1e-12) continue;
        const double first = (grid.origin[axis] - grid.spacing[axis] * 0.5 - source.coordinate[axis]) / delta;
        const double last = (grid.origin[axis] + (grid.dimensions[axis] - 0.5) * grid.spacing[axis] - source.coordinate[axis]) / delta;
        alphaMin = fmax(alphaMin, fmin(first, last));
        alphaMax = fmin(alphaMax, fmax(first, last));
        if (alphaMin >= alphaMax) return;
        const double low = source.coordinate[axis] + (delta > 0.0 ? alphaMin : alphaMax) * delta;
        const double high = source.coordinate[axis] + (delta > 0.0 ? alphaMax : alphaMin) * delta;
        const int begin = max(0, static_cast<int>(ceil((low - grid.origin[axis] + grid.spacing[axis] * 0.5) / grid.spacing[axis])));
        const int end = min(grid.dimensions[axis], static_cast<int>(floor((high - grid.origin[axis] + grid.spacing[axis] * 0.5) / grid.spacing[axis])));
        increment[axis] = delta > 0.0 ? 1 : -1;
        currentPlane[axis] = delta > 0.0 ? begin : end;
        lastPlane[axis] = delta > 0.0 ? end : begin;
        if (begin <= end) hasPlanes = true;
        else { currentPlane[axis] = 0; lastPlane[axis] = -1; increment[axis] = 1; }
    }
    if (!hasPlanes) return;
    bool takeMin = true, takeMax = true, hasPrevious = false;
    double previous = 0.0, cumulative = 0.0, midpointDepth = 0.0;
    const int maxSteps = grid.dimensions[0] + grid.dimensions[1] + grid.dimensions[2] + 8;
    for (int step = 0; step < maxSteps; ++step) {
        double next = CUDART_INF;
        int selected = -1;
        for (int axis = 0; axis < 3; ++axis) {
            if ((increment[axis] > 0 && currentPlane[axis] > lastPlane[axis]) ||
                (increment[axis] < 0 && currentPlane[axis] < lastPlane[axis])) continue;
            const double crossing = (grid.origin[axis] + (currentPlane[axis] - 0.5) * grid.spacing[axis] -
                source.coordinate[axis]) / direction[axis];
            if (crossing < next) { next = crossing; selected = axis; }
        }
        if (takeMin && alphaMin < next) { next = alphaMin; selected = 3; }
        if (takeMax && alphaMax < next) { next = alphaMax; selected = 4; }
        if (selected < 0) break;
        if (selected < 3) currentPlane[selected] += increment[selected];
        else if (selected == 3) takeMin = false;
        else takeMax = false;
        if (hasPrevious && fabs(next - previous) < 1e-12) continue;
        if (hasPrevious) {
            const double length = (next - previous) * distance;
            const double alpha = 0.5 * (next + previous);
            double world[3];
            for (int axis = 0; axis < 3; ++axis) world[axis] = source.coordinate[axis] + alpha * direction[axis];
            int voxel[3];
            bool inside = length >= 1e-12;
            for (int axis = 0; axis < 3; ++axis) {
                double coordinate = 0.0;
                for (int component = 0; component < 3; ++component)
                    coordinate += grid.inverse[axis * 3 + component] * (world[component] - grid.origin[component]);
                coordinate /= grid.spacing[axis];
                voxel[axis] = static_cast<int>(round(coordinate));
                inside = inside && voxel[axis] >= 0 && voxel[axis] < grid.dimensions[axis];
            }
            if (inside) {
                const size_t index = static_cast<size_t>(voxel[0]) * grid.dimensions[1] + voxel[1] +
                    static_cast<size_t>(voxel[2]) * grid.dimensions[0] * grid.dimensions[1];
                const double contribution = length * density[index];
                midpointDepth = cumulative + 0.5 * contribution;
                cumulative += contribution;
            }
        }
        previous = next;
        hasPrevious = true;
    }
    output[ray] = midpointDepth;
}

} // namespace

struct CudaRayTracer::Impl {
    TraceGrid grid;
    int device = 0;
    Buffer<double> density;
    ~Impl() { cudaSetDevice(device); }
};

CudaRayTracer::CudaRayTracer(const Grid& grid, const double* density) : m_impl(std::make_unique<Impl>()) {
    check(cudaGetDevice(&m_impl->device));
    const auto dimensions = grid.getDimensions();
    const auto spacing = grid.getSpacing();
    const auto origin = grid.getOrigin();
    const auto row = grid.getRowDirection(), column = grid.getColumnDirection(), slice = grid.getSliceDirection();
    const Mat3 inverseDirection = inverse(Mat3(column[0], row[0], slice[0], column[1], row[1], slice[1],
        column[2], row[2], slice[2]));
    for (int axis = 0; axis < 3; ++axis) {
        m_impl->grid.dimensions[axis] = static_cast<int>(dimensions[axis]);
        m_impl->grid.spacing[axis] = spacing[axis];
        m_impl->grid.origin[axis] = origin[axis];
    }
    for (int component = 0; component < 3; ++component) {
        Vec3 unit{0.0, 0.0, 0.0};
        unit[component] = 1.0;
        const auto transformed = inverseDirection * unit;
        for (int axis = 0; axis < 3; ++axis) m_impl->grid.inverse[axis * 3 + component] = transformed[axis];
    }
    m_impl->density.allocate(grid.getNumVoxels());
    check(cudaMemcpy(m_impl->density.data, density, grid.getNumVoxels() * sizeof(double), cudaMemcpyHostToDevice));
}

CudaRayTracer::~CudaRayTracer() = default;

std::vector<double> CudaRayTracer::depths(const Vec3& source, const std::vector<Vec3>& targets) {
    check(cudaSetDevice(m_impl->device));
    std::vector<double> result(targets.size());
    if (targets.empty()) return result;
    Point origin{{source[0], source[1], source[2]}};
    std::vector<Point> points;
    points.reserve(targets.size());
    for (const auto& target : targets) points.push_back({{target[0], target[1], target[2]}});
    Buffer<Point> deviceTargets;
    Buffer<double> deviceDepths;
    deviceTargets.allocate(points.size());
    deviceDepths.allocate(points.size());
    check(cudaMemcpy(deviceTargets.data, points.data(), points.size() * sizeof(Point), cudaMemcpyHostToDevice));
    traceDepths<<<static_cast<unsigned>((points.size() + 127) / 128), 128>>>(m_impl->grid, origin,
        deviceTargets.data, points.size(), m_impl->density.data, deviceDepths.data);
    check(cudaGetLastError());
    check(cudaMemcpy(result.data(), deviceDepths.data, result.size() * sizeof(double), cudaMemcpyDeviceToHost));
    return result;
}

} // namespace optirad