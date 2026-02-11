default_target: all
all: examples test

# Compiler flag for OpenMP
OPENMP = -fopenmp

# Disable OpenMP if DISABLE_OPENMP is defined
ifdef DISABLE_OPENMP
	OPENMP = -DTHEORETICA_DISABLE_OPENMP
endif


# Compiler flags
CXXFLAGS = -std=c++14 -I./src/ -Wall ${OPENMP}
CHEBYSHEV_SRC = ./test/chebyshev/src


# Test programs
TEST_SOURCES := $(wildcard test/prec/test_*.cpp)
TEST_BINARIES := $(TEST_SOURCES:.cpp=)

test_%: test/prec/test_%.cpp
	@echo Compiling $(subst test_,,$@) module test unit...
	@g++ $< ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/prec/$@
	@./test/prec/$@

test: $(subst test/prec/,,$(TEST_BINARIES))
	@echo All tests have been compiled and executed.


# Example programs
EXAMPLE_SOURCES := $(wildcard examples/*.cpp)
EXAMPLE_BINARIES := $(EXAMPLE_SOURCES:.cpp=)


# Compile a simple example program
example:
	@echo Compiling quickstart example program...
	@g++ examples/example.cpp ${CXXFLAGS} -o ./examples/example

example_%: examples/%.cpp
	@echo Compiling $(subst example_,,$@) example program...
	@g++ $< ${CXXFLAGS} -o ./examples/$(subst example_,,$@)

examples: example $(subst examples/,example_,$(EXAMPLE_BINARIES))
	@echo All example programs have been compiled.


# Benchmarks
BENCHMARK_SOURCES := $(wildcard test/benchmark/benchmark_*.cpp)
BENCHMARK_BINARIES := $(BENCHMARK_SOURCES:.cpp=)

benchmark_%: test/benchmark/benchmark_%.cpp
	@echo Compiling $(subst benchmark_,,$@) benchmark program...
	@g++ $< ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o test/benchmark/$@
	@./test/benchmark/$@
	@echo All benchmarks have been compiled and executed.

benchmark: $(subst test/benchmark/,,$(BENCHMARK_BINARIES))


# Clean all directories from CSV and EXE files
clean:
	@rm -f ./*.exe
	@rm -f ./test/prec/*.csv
	@rm -f ./test/benchmark/*.csv
	@rm -f ./examples/*.exe
	@rm -f ./test/prec/*.exe
	@rm -f ./test/benchmark/*.exe
