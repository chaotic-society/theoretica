default_target: example
.PHONY: all example test_real_analysis test_algebra autodiff hamiltonian error_propagation stats dist_sample min_grad benchmark clean
all: example test examples

CXXFLAGS = -std=c++14 -I./src/ -Wall -I./include/

example:
	@echo Compiling main example program...
	@g++ examples/example.cpp ${CXXFLAGS} -o ./example


# Tests

test_algebra:
ifndef windows_build
	@echo Compiling linear algebra test cases...
	@g++ test/test_algebra.cpp ${CXXFLAGS} -o test_algebra
	@./test_algebra
endif

test_real_analysis:
ifndef windows_build
	@echo Compiling real analysis test cases...
	@g++ test/test_real_analysis.cpp ${CXXFLAGS} -o test_real_analysis
	@./test_real_analysis
endif

test: test_real_analysis test_algebra


# Example programs

autodiff:
	@echo Compiling \"autodiff\" example...
	@g++ examples/autodiff.cpp ${CXXFLAGS} -o ./autodiff

hamiltonian:
	@echo Compiling \"hamiltonian\" example...
	@g++ examples/hamiltonian.cpp ${CXXFLAGS} -o ./hamiltonian

error_propagation:
	@echo Compiling \"error_propagation\" example...
	@g++ examples/error_propagation.cpp ${CXXFLAGS} -o ./error_propagation

stats:
	@echo Compiling \"stats\" example...
	@g++ examples/stats.cpp ${CXXFLAGS} -o ./stats

dist_sample:
	@echo Compiling \"dist_sample\" example...
	@g++ examples/dist_sample.cpp ${CXXFLAGS} -o ./dist_sample

montecarlo_comparison:
	@echo Compiling \"montecarlo_comparison\" example...
	@g++ examples/montecarlo_comparison.cpp ${CXXFLAGS} -o ./montecarlo_comparison

min_grad:
	@echo Compiling \"min_grad\" example...
	@g++ examples/min_grad.cpp ${CXXFLAGS} -o ./min_grad

examples: autodiff hamiltonian error_propagation stats dist_sample montecarlo_comparison min_grad


# Benchmarks

benchmark_real_analysis:
	@echo Compiling real analysis benchmark...
	@g++ benchmark/benchmark_real_analysis.cpp ${CXXFLAGS} -O0 -o benchmark/benchmark_real_analysis
	@./benchmark/benchmark_real_analysis

benchmark_algebra:
	@echo Compiling algebra benchmark...
	@g++ benchmark/benchmark_algebra.cpp ${CXXFLAGS} -O0 -o benchmark/benchmark_algebra
	@./benchmark/benchmark_algebra

benchmark: benchmark_real_analysis benchmark_algebra
