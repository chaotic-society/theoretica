
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


# Compile a simple example program
example:
	@echo Compiling quickstart example program...
	@g++ examples/example.cpp ${CXXFLAGS} -o ./examples/example


# Tests

test_algebra:
	@echo Compiling linear algebra test cases...
	@g++ test/test_algebra.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_algebra
	@./test/test_algebra


test_core:
	@echo Compiling core test cases...
	@g++ test/test_core.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_core
	@./test/test_core


test_complex:
	@echo Compiling complex test cases...
	@g++ test/test_complex.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_complex
	@./test/test_complex


test_autodiff:
	@echo Compiling autodiff test cases...
	@g++ test/test_autodiff.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_autodiff
	@./test/test_autodiff


test_calculus:
	@echo Compiling calculus test cases...
	@g++ test/test_calculus.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_calculus
	@./test/test_calculus


test_polynomial:
	@echo Compiling polynomial test cases...
	@g++ test/test_polynomial.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_polynomial
	@./test/test_polynomial


test_interpolation:
	@echo Compiling interpolation test cases...
	@g++ test/test_interpolation.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_interpolation
	@./test/test_interpolation


test_optimization:
	@echo Compiling optimization test cases...
	@g++ test/test_optimization.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_optimization
	@./test/test_optimization


test_pseudorandom:
	@echo Compiling pseudorandom test cases...
	@g++ test/test_pseudorandom.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_pseudorandom
	@./test/test_pseudorandom


test_statistics:
	@echo Compiling statistics test cases...
	@g++ test/test_statistics.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_statistics
	@./test/test_statistics


test_template:
	@g++ test/test_template.cpp  ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_template


test_signal:
	@echo Compiling signal test cases...
	@g++ test/test_signal.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -o test/test_signal
	@./test/test_signal


# Compile all test programs and run them
#test: test_core test_algebra test_complex test_autodiff test_calculus test_polynomial test_interpolation test_optimization test_pseudorandom test_statistics
test: test_core test_complex test_statistics test_calculus test_signal

# Example programs

autodiff:
	@echo Compiling \"Automatic differentiation\" example...
	@g++ examples/autodiff.cpp ${CXXFLAGS} -o ./examples/autodiff

hamiltonian:
	@echo Compiling \"Hamiltonian\" example...
	@g++ examples/hamiltonian.cpp ${CXXFLAGS} -o ./examples/hamiltonian

error_propagation:
	@echo Compiling \"Error propagation\" example...
	@g++ examples/error_propagation.cpp ${CXXFLAGS} -o ./examples/error_propagation

statistics:
	@echo Compiling \"Statistics\" example...
	@g++ examples/statistics.cpp ${CXXFLAGS} -o ./examples/statistics

sampling:
	@echo Compiling \"Sampling distributions\" example...
	@g++ examples/sampling.cpp ${CXXFLAGS} -o ./examples/sampling

montecarlo_integral:
	@echo Compiling \"Montecarlo comparison\" example...
	@g++ examples/montecarlo_integral.cpp ${CXXFLAGS} -o ./examples/montecarlo_integral

gradient_descent:
	@echo Compiling \"Multivariate minimization\" example...
	@g++ examples/gradient_descent.cpp ${CXXFLAGS} -o ./examples/gradient_descent

logfit:
	@echo Compiling \"Log fit\" example...
	@g++ examples/logfit.cpp ${CXXFLAGS} -o ./examples/logfit

attractor:
	@echo Compiling \"Attractor\" example...
	@g++ examples/attractor.cpp ${CXXFLAGS} -o ./examples/attractor

histogram:
	@echo Compiling \"Histogram\" example...
	@g++ examples/histogram.cpp ${CXXFLAGS} -o ./examples/histogram

random_walk:
	@echo Compiling \"Random walk\" example...
	@g++ examples/random_walk.cpp ${CXXFLAGS} -o ./examples/random_walk

# Compile all examples
examples: example autodiff hamiltonian error_propagation \
		  statistics sampling montecarlo_integral gradient_descent \
		  logfit attractor histogram random_walk


# Benchmarks

benchmark_real_analysis:
	@echo Compiling real analysis benchmark...
	@g++ benchmark/benchmark_real_analysis.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_real_analysis
	@./benchmark/benchmark_real_analysis

benchmark_algebra:
	@echo Compiling algebra benchmark...
	@g++ benchmark/benchmark_algebra.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_algebra
	@./benchmark/benchmark_algebra

benchmark_dataset:
	@echo Compiling dataset benchmark...
	@g++ benchmark/benchmark_dataset.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_dataset
	@./benchmark/benchmark_dataset

benchmark_pseudorandom:
	@echo Compiling pseudorandom benchmark...
	@g++ benchmark/benchmark_pseudorandom.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_pseudorandom
	@./benchmark/benchmark_pseudorandom

benchmark_vectorized:
	@echo Compiling vectorized benchmark...
	@g++ benchmark/benchmark_vectorized.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_vectorized
	@./benchmark/benchmark_vectorized

benchmark_template:
	@g++ benchmark/benchmark_template.cpp ${CXXFLAGS} -I${CHEBYSHEV_SRC} -O0 -o benchmark/benchmark_template

# Compile all benchmarks and run them
#benchmark: benchmark_real_analysis benchmark_algebra benchmark_dataset benchmark_pseudorandom benchmark_vectorized
benchmark: benchmark_real_analysis benchmark_dataset

# Clean all directories from CSV and EXE files
clean:
	@rm -f ./test/*.csv
	@rm -f ./benchmark/*.csv
	@rm -f ./*.exe
	@rm -f ./examples/*.exe
	@rm -f ./test/*.exe
	@rm -f ./benchmark/*.exe
