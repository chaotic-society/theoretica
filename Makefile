default_target: example
.PHONY: all example test autodiff autodiff_hamiltonian
all: example test

CXXFLAGS = -O2 -DUROBORO_INCLUDE_ALL -DUROBORO_X86 -I./src/

example:
	g++ examples/example.cpp ${CXXFLAGS} -o ./example

test:
	g++ test/test_real_analysis.cpp ${CXXFLAGS} -o test/test_real_analysis
	g++ test/test_algebra.cpp ${CXXFLAGS} -o test/test_algebra
	./test/test_real_analysis
	./test/test_algebra

autodiff:
	g++ examples/autodiff.cpp ${CXXFLAGS} -o ./autodiff

autodiff_hamiltonian:
	g++ examples/autodiff_hamiltonian.cpp ${CXXFLAGS} -o ./autodiff_hamiltonian
