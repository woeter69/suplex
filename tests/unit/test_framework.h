#pragma once

#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <functional>

namespace suplex::test {

struct TestCase {
    std::string name;
    std::function<void()> func;
};

inline std::vector<TestCase>& get_registry() {
    static std::vector<TestCase> registry;
    return registry;
}

inline bool register_test(const std::string& name, std::function<void()> func) {
    get_registry().push_back({name, func});
    return true;
}

#define REGISTER_TEST(test_func) \
    static bool _reg_##test_func = ::suplex::test::register_test(#test_func, test_func)

#define TEST_ASSERT(cond) \
    do { \
        if (!(cond)) { \
            std::cerr << "Assertion failed: " #cond << " at " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::exit(1); \
        } \
    } while (0)

#define TEST_ASSERT_NEAR(val1, val2, eps) \
    do { \
        double _diff = std::abs((val1) - (val2)); \
        if (_diff > (eps)) { \
            std::cerr << "Assertion failed: |" #val1 " - " #val2 "| = " << _diff << " > " << eps \
                      << " (" << (val1) << " vs " << (val2) << ") at " << __FILE__ << ":" << __LINE__ << std::endl; \
            std::exit(1); \
        } \
    } while (0)

} // namespace suplex::test
