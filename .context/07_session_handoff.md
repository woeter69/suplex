# Active Session Snapshot

- **Branch:** `p4`
- **Role completed:** Person 4 — I/O, APIs, GPU integration, CLI, build, and benchmarks.
- **Implemented:** fixed/free MPS and LP readers, solution writer, C++ facade, C ABI,
  pybind11 bindings, CLI, CUDA/CPU SpMV and PCG paths, benchmark runner/reference
  scaffolding, dataset download scripts, sample fixtures, and root CMake integration.
- **Integration fixes:** added the optional `LPSolver::get_farkas_ray()` contract required
  by P3 and made the previously unwired P3 sources clean under `-Werror`.
- **Verification:** 73/73 unit and integration tests pass; a pure C11 client compiles,
  links, and solves successfully; the complete non-optional
  C++ source tree, CLI, and benchmark runner compile with GCC under
  `-std=c++20 -Wall -Wextra -Werror`; sanitizer run passes with LeakSanitizer disabled
  because LSan is unavailable under the sandbox tracer.
- **Environment limitation:** CMake, Make, and Ninja are not installed in this
  environment, so the equivalent source sets were compiled directly with GCC. CUDA
  and pybind11 are optional and their toolchains are not installed here.
- **Upstream gap:** `origin/p2` has context only and no `src/ipm/` implementation.
  IPM and QP requests therefore return an explicit unavailable diagnostic.
- **Presolve caveat:** P3's postsolve stack lacks a complete original-to-reduced column
  map and misclassifies a basic binary smoke problem. LP and MILP currently solve in
  original space to preserve correct status, primal values, and objective values.

## Next Integration Actions

1. Merge the actual P2 IPM/QP implementation and enable its dispatch in `Suplex`.
2. Complete P3 postsolve column mapping, then enable continuous-LP presolve dispatch.
3. Validate optional builds on hosts with CMake 3.24+, pybind11, and CUDA 12+.
4. Download full benchmark corpora and populate the complete reference CSV files.
