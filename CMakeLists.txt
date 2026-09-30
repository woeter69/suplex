cmake_minimum_required(VERSION 3.24)
project(suplex VERSION 0.1.0 LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

option(SUPLEX_BUILD_TESTS "Build the test suite" ON)
option(SUPLEX_BUILD_BENCHMARKS "Build the benchmark runner" ON)
option(SUPLEX_BUILD_PYTHON "Build Python bindings" OFF)
option(SUPLEX_ENABLE_GPU "Enable CUDA acceleration" OFF)
option(SUPLEX_BUILD_P3_GTESTS "Build Person 3 GoogleTest suites when available" ON)

if(CMAKE_CXX_COMPILER_ID MATCHES "Clang|GNU")
    add_compile_options(-Wall -Wextra -Werror)
endif()

add_library(suplex_core src/core/sparse_matrix.cpp src/core/problem.cpp src/core/solution.cpp)
target_include_directories(suplex_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(suplex_core PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(suplex_simplex
    src/simplex/basis.cpp src/simplex/lu_factor.cpp src/simplex/lu_update.cpp
    src/simplex/pricing.cpp src/simplex/primal_simplex.cpp src/simplex/dual_simplex.cpp)
target_link_libraries(suplex_simplex PUBLIC suplex_core)
target_include_directories(suplex_simplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(suplex_simplex PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(suplex_milp
    src/milp/presolve.cpp src/milp/node.cpp src/milp/branching.cpp src/milp/cuts.cpp
    src/milp/heuristics.cpp src/milp/branch_bound.cpp src/milp/conflict.cpp src/milp/milp_solver.cpp)
target_link_libraries(suplex_milp PUBLIC suplex_simplex suplex_core)
target_include_directories(suplex_milp PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(suplex_milp PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(suplex_io src/io/mps_parser.cpp src/io/lp_parser.cpp src/io/sol_writer.cpp)
target_link_libraries(suplex_io PUBLIC suplex_core)
target_include_directories(suplex_io PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(suplex_io PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_library(suplex_gpu src/gpu/gpu_manager.cpp)
target_link_libraries(suplex_gpu PUBLIC suplex_core)
target_include_directories(suplex_gpu PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(suplex_gpu PROPERTIES POSITION_INDEPENDENT_CODE ON)
if(SUPLEX_ENABLE_GPU)
    enable_language(CUDA)
    find_package(CUDAToolkit REQUIRED)
    target_sources(suplex_gpu PRIVATE src/gpu/cuda_spmv.cu src/gpu/cuda_pcg.cu)
    target_compile_definitions(suplex_gpu PUBLIC SUPLEX_GPU_ENABLED)
    target_link_libraries(suplex_gpu PRIVATE CUDA::cudart)
    set_target_properties(suplex_gpu PROPERTIES CUDA_STANDARD 17 CUDA_STANDARD_REQUIRED ON)
endif()

add_library(suplex SHARED src/api/suplex_solver.cpp src/api/suplex_api.cpp)
target_link_libraries(suplex PUBLIC suplex_io suplex_milp suplex_simplex suplex_gpu suplex_core)
target_include_directories(suplex PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/src/api)
set_target_properties(suplex PROPERTIES VERSION ${PROJECT_VERSION} SOVERSION 0)

add_executable(suplex_cli cli/main.cpp)
target_link_libraries(suplex_cli PRIVATE suplex)
set_target_properties(suplex_cli PROPERTIES OUTPUT_NAME suplex)

if(SUPLEX_BUILD_TESTS)
    enable_testing()
    add_executable(suplex_test_runner
        tests/unit/test_runner.cpp tests/unit/test_sparse_matrix.cpp tests/unit/test_lu.cpp
        tests/unit/test_simplex.cpp tests/unit/test_boundary.cpp tests/integration/test_integration.cpp)
    target_link_libraries(suplex_test_runner PRIVATE suplex)
    target_compile_definitions(suplex_test_runner PRIVATE SUPLEX_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
    add_test(NAME suplex_unit_and_integration_tests COMMAND suplex_test_runner)
    add_test(NAME suplex_cli_help COMMAND suplex_cli --help)
    add_test(NAME suplex_cli_solve COMMAND suplex_cli --algorithm primal --solution
        ${CMAKE_CURRENT_BINARY_DIR}/tiny.sol ${CMAKE_CURRENT_SOURCE_DIR}/data/examples/tiny.mps)
    set_tests_properties(suplex_cli_solve PROPERTIES PASS_REGULAR_EXPRESSION "Status: OPTIMAL")
    add_test(NAME suplex_cli_rejects_unknown_option COMMAND suplex_cli --not-an-option value)
    set_tests_properties(suplex_cli_rejects_unknown_option PROPERTIES WILL_FAIL TRUE)
    add_executable(suplex_c_api_test tests/integration/test_c_api.c)
    target_link_libraries(suplex_c_api_test PRIVATE suplex)
    add_test(NAME suplex_c_api_client COMMAND suplex_c_api_test)
    add_test(NAME suplex_download_script_contracts
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_scripts.sh
        ${CMAKE_CURRENT_SOURCE_DIR})
    add_test(NAME suplex_cli_full_contract
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_cli.sh
        $<TARGET_FILE:suplex_cli> ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_BINARY_DIR}/cli-test-output)

    if(SUPLEX_BUILD_P3_GTESTS)
        find_package(GTest QUIET)
        if(GTest_FOUND)
            foreach(p3_suite IN ITEMS presolve branching cuts heuristics)
                add_executable(suplex_p3_${p3_suite}_test tests/unit/milp/test_${p3_suite}.cpp)
                target_link_libraries(suplex_p3_${p3_suite}_test PRIVATE suplex GTest::gtest)
                target_include_directories(suplex_p3_${p3_suite}_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
                add_test(NAME suplex_p3_${p3_suite} COMMAND suplex_p3_${p3_suite}_test)
            endforeach()
        else()
            message(STATUS "GoogleTest not found: Person 3 standalone suites are source-only")
        endif()
    endif()
endif()

if(SUPLEX_BUILD_BENCHMARKS)
    add_executable(suplex_benchmark tests/benchmarks/benchmark_runner.cpp)
    target_link_libraries(suplex_benchmark PRIVATE suplex)
    if(SUPLEX_BUILD_TESTS)
        add_test(NAME suplex_benchmark_smoke
            COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_benchmark.sh
            $<TARGET_FILE:suplex_benchmark> ${CMAKE_CURRENT_SOURCE_DIR}
            ${CMAKE_CURRENT_BINARY_DIR}/benchmark-test-output)
    endif()
endif()

if(SUPLEX_BUILD_PYTHON)
    find_package(pybind11 CONFIG REQUIRED)
    pybind11_add_module(pysuplex src/api/python/bindings.cpp)
    target_link_libraries(pysuplex PRIVATE suplex)
    if(SUPLEX_BUILD_TESTS)
        find_package(Python3 COMPONENTS Interpreter REQUIRED)
        add_test(NAME suplex_python_api
            COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/tests/integration/test_python_api.py
            ${CMAKE_CURRENT_SOURCE_DIR})
        set_tests_properties(suplex_python_api PROPERTIES
            ENVIRONMENT "PYTHONPATH=$<TARGET_FILE_DIR:pysuplex>")
    endif()
endif()

install(TARGETS suplex suplex_cli)
install(FILES src/api/suplex.h DESTINATION include)
