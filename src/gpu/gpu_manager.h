#pragma once

#include "src/core/sparse_matrix.h"

#include <cstddef>
#include <string>

namespace suplex {

/** Singleton CUDA capability manager with transparent CPU fallbacks. */
class GPUManager {
public:
    static GPUManager& instance();
    bool initialize(int device_id = 0);
    bool is_available() const { return available_; }
    int device_id() const { return device_id_; }
    const std::string& last_error() const { return error_; }
    void synchronize();
    void* device_malloc(std::size_t bytes);
    void device_free(void* pointer);
    bool copy_to_device(void* destination, const void* source, std::size_t bytes);
    bool copy_to_host(void* destination, const void* source, std::size_t bytes);
    void spmv_csr(const SparseMatrixCSR& matrix, const Real* x, Real* y);
    bool pcg_solve(const SparseMatrixCSC& matrix, const Real* b, Real* x,
                   Real tolerance, Index max_iterations);

private:
    GPUManager() = default;
    bool available_ = false;
    int device_id_ = -1;
    std::string error_;
};

} // namespace suplex
