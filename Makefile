CXX ?= g++
CC ?= gcc
AR ?= ar
CPPFLAGS := -I.
CXXFLAGS ?= -O2 -g
CXXFLAGS += -std=c++20 -Wall -Wextra -Werror -fPIC
CFLAGS ?= -O2 -g
CFLAGS += -std=c11 -Wall -Wextra -Werror

BUILD_DIR := build_make
LIB_DIR := $(BUILD_DIR)/lib
BIN_DIR := $(BUILD_DIR)/bin

LIB_SOURCES := $(wildcard src/core/*.cpp) \
               $(wildcard src/simplex/*.cpp) \
               $(wildcard src/milp/*.cpp) \
               $(wildcard src/io/*.cpp) \
               $(filter-out src/api/python/bindings.cpp,$(wildcard src/api/*.cpp)) \
               src/gpu/gpu_manager.cpp
LIB_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(LIB_SOURCES))
LIBRARY := $(LIB_DIR)/libsuplex.a

TEST_SOURCES := tests/unit/test_runner.cpp \
                tests/unit/test_sparse_matrix.cpp \
                tests/unit/test_lu.cpp \
                tests/unit/test_simplex.cpp \
                tests/unit/test_boundary.cpp \
                tests/integration/test_integration.cpp
TEST_OBJECTS := $(patsubst %.cpp,$(BUILD_DIR)/%.o,$(TEST_SOURCES))

.PHONY: all test benchmark clean

all: $(BIN_DIR)/suplex

$(LIBRARY): $(LIB_OBJECTS)
	@mkdir -p $(dir $@)
	$(AR) rcs $@ $^

$(BIN_DIR)/suplex: $(BUILD_DIR)/cli/main.o $(LIBRARY)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/suplex_tests: CPPFLAGS += -DSUPLEX_SOURCE_DIR=\"$(CURDIR)\"
$(BIN_DIR)/suplex_tests: $(TEST_OBJECTS) $(LIBRARY)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@

$(BIN_DIR)/suplex_c_api_test: $(BUILD_DIR)/tests/integration/test_c_api.o $(LIBRARY)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -lm -o $@

$(BIN_DIR)/suplex_benchmark: $(BUILD_DIR)/tests/benchmarks/benchmark_runner.o $(LIBRARY)
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $^ -o $@

test: $(BIN_DIR)/suplex_tests $(BIN_DIR)/suplex_c_api_test
	$(BIN_DIR)/suplex_tests
	$(BIN_DIR)/suplex_c_api_test
	$(BIN_DIR)/suplex --help >/dev/null
	$(BIN_DIR)/suplex --algorithm primal data/examples/tiny.mps | grep -q "Status: OPTIMAL"

benchmark: $(BIN_DIR)/suplex_benchmark

$(BUILD_DIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

$(BUILD_DIR)/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(BUILD_DIR)
