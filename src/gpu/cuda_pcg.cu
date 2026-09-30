#include "gpu_kernels.h"

#include <cuda_runtime.h>

#include <cmath>
#include <vector>

namespace suplex {
namespace {

__global__ void pcg_spmv(int rows, const int* row_start, const int* columns,
                         const double* values, const double* x, double* y) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= rows) return;
    double sum = 0.0;
    for (int k = row_start[row]; k < row_start[row + 1]; ++k) sum += values[k] * x[columns[k]];
    y[row] = sum;
}

__global__ void initialize_residual(int n, const double* b, const double* ax,
                                    const double* diagonal, double* r, double* z, double* p) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) { r[i] = b[i] - ax[i]; z[i] = r[i] / diagonal[i]; p[i] = z[i]; }
}

__global__ void dot_kernel(int n, const double* a, const double* b, double* result) {
    __shared__ double partial[256];
    const int tid = threadIdx.x;
    const int i = blockIdx.x * blockDim.x + tid;
    partial[tid] = i < n ? a[i] * b[i] : 0.0;
    __syncthreads();
    for (int stride = blockDim.x / 2; stride > 0; stride >>= 1) {
        if (tid < stride) partial[tid] += partial[tid + stride];
        __syncthreads();
    }
    if (tid == 0) atomicAdd(result, partial[0]);
}

__global__ void update_x_r(int n, double alpha, const double* p, const double* ap,
                           double* x, double* r) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) { x[i] += alpha * p[i]; r[i] -= alpha * ap[i]; }
}

__global__ void apply_jacobi(int n, const double* diagonal, const double* r, double* z) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) z[i] = r[i] / diagonal[i];
}

__global__ void update_direction(int n, double beta, const double* z, double* p) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) p[i] = z[i] + beta * p[i];
}

bool device_dot(int n, const double* a, const double* b, double* scalar, double& host) {
    if (cudaMemset(scalar, 0, sizeof(double)) != cudaSuccess) return false;
    constexpr int threads = 256;
    dot_kernel<<<(n + threads - 1) / threads, threads>>>(n, a, b, scalar);
    return cudaGetLastError() == cudaSuccess &&
           cudaMemcpy(&host, scalar, sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
}

} // namespace

bool cuda_pcg_host(const SparseMatrixCSC& matrix, const Real* b, Real* x,
                   Real tolerance, Index max_iterations) {
    if (matrix.num_rows != matrix.num_cols || matrix.num_rows <= 0) return false;
    const SparseMatrixCSR csr = matrix.to_csr();
    const int n = csr.num_rows;
    std::vector<double> diagonal(static_cast<std::size_t>(n), 1.0);
    for (int i = 0; i < n; ++i) {
        const double value = csr.get(i, i);
        if (std::abs(value) > EPS_ZERO) diagonal[static_cast<std::size_t>(i)] = value;
    }

    int *d_row = nullptr, *d_col = nullptr;
    double *d_value = nullptr, *d_b = nullptr, *d_x = nullptr, *d_r = nullptr;
    double *d_z = nullptr, *d_p = nullptr, *d_ap = nullptr, *d_diag = nullptr, *d_scalar = nullptr;
    auto allocate = [](void** pointer, std::size_t bytes) { return cudaMalloc(pointer, bytes) == cudaSuccess; };
    const std::size_t vector_bytes = static_cast<std::size_t>(n) * sizeof(double);
    bool ok = allocate(reinterpret_cast<void**>(&d_row), csr.row_start.size() * sizeof(int)) &&
        allocate(reinterpret_cast<void**>(&d_col), csr.col_index.size() * sizeof(int)) &&
        allocate(reinterpret_cast<void**>(&d_value), csr.values.size() * sizeof(double)) &&
        allocate(reinterpret_cast<void**>(&d_b), vector_bytes) && allocate(reinterpret_cast<void**>(&d_x), vector_bytes) &&
        allocate(reinterpret_cast<void**>(&d_r), vector_bytes) && allocate(reinterpret_cast<void**>(&d_z), vector_bytes) &&
        allocate(reinterpret_cast<void**>(&d_p), vector_bytes) && allocate(reinterpret_cast<void**>(&d_ap), vector_bytes) &&
        allocate(reinterpret_cast<void**>(&d_diag), vector_bytes) && allocate(reinterpret_cast<void**>(&d_scalar), sizeof(double));
    if (ok) {
        ok = cudaMemcpy(d_row, csr.row_start.data(), csr.row_start.size() * sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess &&
            cudaMemcpy(d_col, csr.col_index.data(), csr.col_index.size() * sizeof(int), cudaMemcpyHostToDevice) == cudaSuccess &&
            cudaMemcpy(d_value, csr.values.data(), csr.values.size() * sizeof(double), cudaMemcpyHostToDevice) == cudaSuccess &&
            cudaMemcpy(d_b, b, vector_bytes, cudaMemcpyHostToDevice) == cudaSuccess &&
            cudaMemcpy(d_x, x, vector_bytes, cudaMemcpyHostToDevice) == cudaSuccess &&
            cudaMemcpy(d_diag, diagonal.data(), vector_bytes, cudaMemcpyHostToDevice) == cudaSuccess;
    }
    constexpr int threads = 256;
    const int blocks = (n + threads - 1) / threads;
    if (ok) {
        pcg_spmv<<<blocks, threads>>>(n, d_row, d_col, d_value, d_x, d_ap);
        initialize_residual<<<blocks, threads>>>(n, d_b, d_ap, d_diag, d_r, d_z, d_p);
        ok = cudaGetLastError() == cudaSuccess;
    }
    double rz = 0.0;
    if (ok) ok = device_dot(n, d_r, d_z, d_scalar, rz);
    bool converged = false;
    for (Index iteration = 0; ok && iteration < max_iterations; ++iteration) {
        pcg_spmv<<<blocks, threads>>>(n, d_row, d_col, d_value, d_p, d_ap);
        double denominator = 0.0;
        ok = device_dot(n, d_p, d_ap, d_scalar, denominator) && std::abs(denominator) > EPS_ZERO;
        if (!ok) break;
        update_x_r<<<blocks, threads>>>(n, rz / denominator, d_p, d_ap, d_x, d_r);
        double norm2 = 0.0;
        ok = device_dot(n, d_r, d_r, d_scalar, norm2);
        if (ok && std::sqrt(norm2) <= tolerance) { converged = true; break; }
        apply_jacobi<<<blocks, threads>>>(n, d_diag, d_r, d_z);
        double next_rz = 0.0;
        ok = device_dot(n, d_r, d_z, d_scalar, next_rz) && std::abs(rz) > EPS_ZERO;
        if (ok) update_direction<<<blocks, threads>>>(n, next_rz / rz, d_z, d_p);
        rz = next_rz;
    }
    if (ok && converged) ok = cudaMemcpy(x, d_x, vector_bytes, cudaMemcpyDeviceToHost) == cudaSuccess;
    cudaFree(d_row); cudaFree(d_col); cudaFree(d_value); cudaFree(d_b); cudaFree(d_x);
    cudaFree(d_r); cudaFree(d_z); cudaFree(d_p); cudaFree(d_ap); cudaFree(d_diag); cudaFree(d_scalar);
    return ok && converged;
}

} // namespace suplex
