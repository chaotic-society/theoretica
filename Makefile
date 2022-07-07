default_target: example
.PHONY: all example test autodiff hamiltonian error_propagation stats dist_sample benchmark
all: example test examples

CXXFLAGS = -std=c++14 -I./src/ -Wall
EXE_FORMAT = 

example:
	@echo Compiling main example program...
	@g++ examples/example.cpp ${CXXFLAGS} -o ./example$(EXE_FORMAT)

test_algebra:
	@echo Compiling linear algebra test cases...
	@g++ test/test_algebra.cpp ${CXXFLAGS} -o test/test_algebra$(EXE_FORMAT)
	@./test/test_algebra$(EXE_FORMAT)

test_real_analysis:
	@echo Compiling real analysis test cases...
	@g++ test/test_real_analysis.cpp ${CXXFLAGS} -o test/test_real_analysis$(EXE_FORMAT)
	@./test/test_real_analysis$(EXE_FORMAT)

test: test_real_analysis test_algebra

autodiff:
	@echo Compiling \"autodiff\" example...
	@g++ examples/autodiff.cpp ${CXXFLAGS} -o ./autodiff$(EXE_FORMAT)

hamiltonian:
	@echo Compiling \"hamiltonian\" example...
	@g++ examples/hamiltonian.cpp ${CXXFLAGS} -o ./hamiltonian$(EXE_FORMAT)

error_propagation:
	@echo Compiling \"error_propagation\" example...
	@g++ examples/error_propagation.cpp ${CXXFLAGS} -o ./error_propagation$(EXE_FORMAT)

stats:
	@echo Compiling \"stats\" example...
	@g++ examples/stats.cpp ${CXXFLAGS} -o ./stats$(EXE_FORMAT)

dist_sample:
	@echo Compiling \"dist_sample\" example...
	@g++ examples/dist_sample.cpp ${CXXFLAGS} -o ./dist_sample$(EXE_FORMAT)

examples: autodiff hamiltonian error_propagation stats dist_sample

benchmark_real_analysis:
	@echo Compiling real analysis benchmark...
	@g++ benchmark/benchmark_real_analysis.cpp ${CXXFLAGS} -O0 -o benchmark/benchmark_real_analysis$(EXE_FORMAT)
	@./benchmark/benchmark_real_analysis$(EXE_FORMAT)

benchmark: benchmark_real_analysis
