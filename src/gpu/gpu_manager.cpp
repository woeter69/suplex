#include "gpu_manager.h"

#include "gpu_kernels.h"

#include <algorithm>
#include <cmath>
#include <vector>

#ifdef SUPLEX_GPU_ENABLED
#include <cuda_runtime_api.h>
#endif

namespace suplex {

GPUManager& GPUManager::instance() {
    static GPUManager manager;
    return manager;
}

bool GPUManager::initialize(int device_id) {
#ifdef SUPLEX_GPU_ENABLED
    int count = 0;
    if (cudaGetDeviceCount(&count) != cudaSuccess || count == 0) {
        error_ = "no CUDA device is available"; available_ = false; return false;
    }
    if (device_id < 0 || device_id >= count || cudaSetDevice(device_id) != cudaSuccess) {
        error_ = "invalid or inaccessible CUDA device"; available_ = false; return false;
    }
    device_id_ = device_id; available_ = true; error_.clear(); return true;
#else
    (void)device_id;
    error_ = "Suplex was built without CUDA support";
    available_ = false;
    return false;
#endif
}

void GPUManager::synchronize() {
#ifdef SUPLEX_GPU_ENABLED
    if (available_ && cudaDeviceSynchronize() != cudaSuccess) error_ = "CUDA synchronization failed";
#endif
}

void* GPUManager::device_malloc(std::size_t bytes) {
#ifdef SUPLEX_GPU_ENABLED
    if (!available_) return nullptr;
    void* pointer = nullptr;
    if (cudaMalloc(&pointer, bytes) != cudaSuccess) { error_ = "CUDA allocation failed"; return nullptr; }
    return pointer;
#else
    (void)bytes; return nullptr;
#endif
}

void GPUManager::device_free(void* pointer) {
#ifdef SUPLEX_GPU_ENABLED
    if (pointer) cudaFree(pointer);
#else
    (void)pointer;
#endif
}

bool GPUManager::copy_to_device(void* destination, const void* source, std::size_t bytes) {
#ifdef SUPLEX_GPU_ENABLED
    return available_ && cudaMemcpy(destination, source, bytes, cudaMemcpyHostToDevice) == cudaSuccess;
#else
    (void)destination; (void)source; (void)bytes; return false;
#endif
}

bool GPUManager::copy_to_host(void* destination, const void* source, std::size_t bytes) {
#ifdef SUPLEX_GPU_ENABLED
    return available_ && cudaMemcpy(destination, source, bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
#else
    (void)destination; (void)source; (void)bytes; return false;
#endif
}

void GPUManager::spmv_csr(const SparseMatrixCSR& matrix, const Real* x, Real* y) {
#ifdef SUPLEX_GPU_ENABLED
    if (available_ && matrix.nnz > 100000 && cuda_spmv_host(matrix, x, y)) return;
#endif
    matrix.multiply_vec(x, y);
}

bool GPUManager::pcg_solve(const SparseMatrixCSC& matrix, const Real* b, Real* x,
                           Real tolerance, Index max_iterations) {
#ifdef SUPLEX_GPU_ENABLED
    if (available_ && matrix.nnz > 100000 && cuda_pcg_host(matrix, b, x, tolerance, max_iterations)) return true;
#endif
    if (matrix.num_rows != matrix.num_cols || tolerance <= 0.0 || max_iterations <= 0) return false;
    const Index n = matrix.num_rows;
    const SparseMatrixCSR csr = matrix.to_csr();
    std::vector<Real> r(static_cast<std::size_t>(n)), p(static_cast<std::size_t>(n));
    std::vector<Real> ap(static_cast<std::size_t>(n)), diagonal(static_cast<std::size_t>(n), 1.0);
    std::vector<Real> z(static_cast<std::size_t>(n));
    for (Index i = 0; i < n; ++i) {
        const Real value = csr.get(i, i);
        if (std::abs(value) > EPS_ZERO) diagonal[static_cast<std::size_t>(i)] = value;
    }
    csr.multiply_vec(x, ap.data());
    Real rz = 0.0;
    for (Index i = 0; i < n; ++i) {
        const std::size_t k = static_cast<std::size_t>(i);
        r[k] = b[k] - ap[k]; z[k] = r[k] / diagonal[k]; p[k] = z[k]; rz += r[k] * z[k];
    }
    for (Index iteration = 0; iteration < max_iterations; ++iteration) {
        csr.multiply_vec(p.data(), ap.data());
        Real denominator = 0.0;
        for (Index i = 0; i < n; ++i) denominator += p[static_cast<std::size_t>(i)] * ap[static_cast<std::size_t>(i)];
        if (std::abs(denominator) <= EPS_ZERO) return false;
        const Real alpha = rz / denominator;
        Real norm2 = 0.0;
        for (Index i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i);
            x[k] += alpha * p[k]; r[k] -= alpha * ap[k]; norm2 += r[k] * r[k];
        }
        if (std::sqrt(norm2) <= tolerance) return true;
        Real next_rz = 0.0;
        for (Index i = 0; i < n; ++i) {
            const std::size_t k = static_cast<std::size_t>(i); z[k] = r[k] / diagonal[k]; next_rz += r[k] * z[k];
        }
        const Real beta = next_rz / rz;
        for (Index i = 0; i < n; ++i) { const std::size_t k = static_cast<std::size_t>(i); p[k] = z[k] + beta * p[k]; }
        rz = next_rz;
    }
    return false;
}

} // namespace suplex
