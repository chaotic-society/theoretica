default_target: example
.PHONY: all example test autodiff hamiltonian error_propagation stats dist_sample
all: example test examples

CXXFLAGS = -DTHEORETICA_INCLUDE_ALL -std=c++14 -I./src/

example:
	@echo Compiling main example program...
	@g++ examples/example.cpp ${CXXFLAGS} -o ./example

test_algebra:
	@echo Compiling linear algebra test cases...
	@g++ test/test_algebra.cpp ${CXXFLAGS} -o test/test_algebra
	@./test/test_algebra

test_real_analysis:
	@echo Compiling real analysis test cases...
	@g++ test/test_real_analysis.cpp ${CXXFLAGS} -o test/test_real_analysis
	@./test/test_real_analysis

test: test_real_analysis test_algebra

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

examples: autodiff hamiltonian error_propagation stats dist_sample
