# Suplex Person 3 — Makefile (direct g++ build, no cmake required)
CXX      := g++
CXXFLAGS := -std=c++20 -Wall -Wextra -Wpedantic -g -O0 \
            -Isrc \
            -fsanitize=address,undefined
LDFLAGS  := -fsanitize=address,undefined

# ── Sources ────────────────────────────────────────────────────────────────────
CORE_SRC := src/core/core_stubs.cpp

MILP_SRC := src/milp/presolve.cpp \
            src/milp/node.cpp \
            src/milp/branching.cpp \
            src/milp/cuts.cpp \
            src/milp/heuristics.cpp \
            src/milp/branch_bound.cpp \
            src/milp/conflict.cpp \
            src/milp/milp_solver.cpp

ALL_SRC := $(CORE_SRC) $(MILP_SRC)

# ── Build directory ────────────────────────────────────────────────────────────
BUILDDIR := build_make

# ── Object files ──────────────────────────────────────────────────────────────
OBJS := $(patsubst %.cpp,$(BUILDDIR)/%.o,$(ALL_SRC))

# ── Libraries (static archive) ────────────────────────────────────────────────
LIB := $(BUILDDIR)/libsuplex_milp.a

# ── Test executables ──────────────────────────────────────────────────────────
TEST_SRC := tests/unit/milp/test_presolve.cpp \
            tests/unit/milp/test_branching.cpp \
            tests/unit/milp/test_cuts.cpp \
            tests/unit/milp/test_heuristics.cpp

# We use a single-file compile for tests (no GTest in env → run syntax-check only)
TEST_BINS := $(patsubst tests/unit/milp/%.cpp,$(BUILDDIR)/tests/%,$(TEST_SRC))

# ── Default target: compile library ───────────────────────────────────────────
.PHONY: all lib tests clean check

all: lib

lib: $(LIB)

$(LIB): $(OBJS)
	@mkdir -p $(dir $@)
	ar rcs $@ $^
	@echo "✅  Library built: $@"

$(BUILDDIR)/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

# ── Syntax-check all sources without linking ──────────────────────────────────
check:
	@echo "Syntax-checking all sources..."
	@for f in $(ALL_SRC); do \
	    echo "  checking $$f ..."; \
	    $(CXX) $(CXXFLAGS) -fsyntax-only $$f && echo "    OK" || echo "    FAIL $$f"; \
	done
	@echo "Done."

# ── Build test objects (syntax check only — no GTest linkage needed) ──────────
check-tests:
	@echo "Syntax-checking test sources..."
	@for f in $(TEST_SRC); do \
	    echo "  checking $$f ..."; \
	    $(CXX) $(CXXFLAGS) -fsyntax-only $$f && echo "    OK" || echo "    FAIL $$f"; \
	done
	@echo "Done."

clean:
	rm -rf $(BUILDDIR)
	@echo "Cleaned."
