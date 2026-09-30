#include "test_framework.h"
#include <iostream>
#include <iomanip>

int main() {
    auto& registry = suplex::test::get_registry();
    std::cout << "========================================" << std::endl;
    std::cout << "Running Suplex Test Suite" << std::endl;
    std::cout << "Total tests registered: " << registry.size() << std::endl;
    std::cout << "========================================" << std::endl;

    int passed = 0;
    int failed = 0;

    for (const auto& test_case : registry) {
        std::cout << "[RUN ] " << test_case.name << " ... " << std::flush;
        try {
            test_case.func();
            std::cout << "[PASS]" << std::endl;
            passed++;
        } catch (const std::exception& ex) {
            std::cout << "[FAIL] (Exception: " << ex.what() << ")" << std::endl;
            failed++;
        } catch (...) {
            std::cout << "[FAIL] (Unknown exception)" << std::endl;
            failed++;
        }
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Test Summary: " << passed << " passed, " << failed << " failed." << std::endl;
    std::cout << "========================================" << std::endl;

    return (failed == 0) ? 0 : 1;
}
