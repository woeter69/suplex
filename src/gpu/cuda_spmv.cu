#include "gpu_kernels.h"

#include <cuda_runtime.h>

namespace suplex {
namespace {

__global__ void spmv_csr_kernel(int rows, const int* row_start, const int* column,
                                const double* values, const double* x, double* y) {
    const int row = blockIdx.x * blockDim.x + threadIdx.x;
    if (row >= rows) return;
    double sum = 0.0;
    for (int k = row_start[row]; k < row_start[row + 1]; ++k) sum += values[k] * x[column[k]];
    y[row] = sum;
}

__global__ void spmv_csr_warp_kernel(int rows, const int* row_start, const int* column,
                                     const double* values, const double* x, double* y) {
    const int lane = threadIdx.x & 31;
    const int row = (blockIdx.x * blockDim.x + threadIdx.x) >> 5;
    if (row >= rows) return;
    double sum = 0.0;
    for (int k = row_start[row] + lane; k < row_start[row + 1]; k += 32) sum += values[k] * x[column[k]];
    for (int offset = 16; offset > 0; offset >>= 1) sum += __shfl_down_sync(0xffffffffU, sum, offset);
    if (lane == 0) y[row] = sum;
}

template <typename T>
bool allocate_copy(T*& device, const T* host, std::size_t count) {
    return cudaMalloc(&device, count * sizeof(T)) == cudaSuccess &&
           cudaMemcpy(device, host, count * sizeof(T), cudaMemcpyHostToDevice) == cudaSuccess;
}

} // namespace

bool cuda_spmv_host(const SparseMatrixCSR& matrix, const Real* x, Real* y) {
    int *row = nullptr, *column = nullptr;
    double *value = nullptr, *device_x = nullptr, *device_y = nullptr;
    const bool allocated = allocate_copy(row, matrix.row_start.data(), matrix.row_start.size()) &&
        allocate_copy(column, matrix.col_index.data(), matrix.col_index.size()) &&
        allocate_copy(value, matrix.values.data(), matrix.values.size()) &&
        allocate_copy(device_x, x, static_cast<std::size_t>(matrix.num_cols)) &&
        cudaMalloc(&device_y, static_cast<std::size_t>(matrix.num_rows) * sizeof(double)) == cudaSuccess;
    if (allocated) {
        constexpr int threads = 256;
        const int blocks = (matrix.num_rows * 32 + threads - 1) / threads;
        spmv_csr_warp_kernel<<<blocks, threads>>>(matrix.num_rows, row, column, value, device_x, device_y);
    }
    const bool ok = allocated && cudaGetLastError() == cudaSuccess &&
        cudaMemcpy(y, device_y, static_cast<std::size_t>(matrix.num_rows) * sizeof(double), cudaMemcpyDeviceToHost) == cudaSuccess;
    cudaFree(row); cudaFree(column); cudaFree(value); cudaFree(device_x); cudaFree(device_y);
    return ok;
}

} // namespace suplex
