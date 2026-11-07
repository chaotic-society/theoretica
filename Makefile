default_target: all
all: examples test


# Configuration

# C++ Standard (c++14, c++17, c++20)
CXX_STANDARD ?= c++20


# OpenMP support
# Define DISABLE_OPENMP to disable OpenMP support.
OPENMP = -fopenmp

ifdef DISABLE_OPENMP
OPENMP = -DTHEORETICA_DISABLE_OPENMP
endif


# HDF5 support
# Change the paths below to match your HDF5 installation.
# If empty, HDF5 support will be disabled.
HDF5_ROOT =
HDF5_INCLUDE = $(HDF5_ROOT)/include/
HDF5_LIB = $(HDF5_ROOT)/lib/
HDF5_FLAGS = 

ifneq ($(strip $(HDF5_ROOT)),)
HDF5_FLAGS = -I"$(HDF5_INCLUDE)" -L"$(HDF5_LIB)" -DTHEORETICA_HAS_HDF5 -lhdf5 -lzlib -lsz
else
$(info HDF5 support is disabled. To enable it, set HDF5_ROOT in the Makefile.)
endif


# Flags

# Compiler flags
CXXFLAGS = -std=${CXX_STANDARD} -I./src/ -Wall ${OPENMP}
CHEBYSHEV_SRC = ./test/chebyshev/src

# Test programs
TEST_SOURCES := $(wildcard test/prec/test_*.cpp)
TEST_BINARIES := $(TEST_SOURCES:.cpp=)


# Targets

test_%: test/prec/test_%.cpp
	@echo Compiling $(subst test_,,$@) module test unit...
	@g++ $< ${CXXFLAGS} -I${CHEBYSHEV_SRC} $(if $(filter test_io,$@),$(HDF5_FLAGS),) -o test/prec/$@
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
	@g++ $< ${CXXFLAGS} $(if $(filter example_hdf5_io,$@),$(HDF5_FLAGS),) -o ./examples/$(subst example_,,$@)

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
