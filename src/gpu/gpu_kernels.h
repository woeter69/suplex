#pragma once

#include "src/core/sparse_matrix.h"

namespace suplex {

#ifdef SUPLEX_GPU_ENABLED
bool cuda_spmv_host(const SparseMatrixCSR& matrix, const Real* x, Real* y);
bool cuda_pcg_host(const SparseMatrixCSC& matrix, const Real* b, Real* x,
                   Real tolerance, Index max_iterations);
#endif

} // namespace suplex
