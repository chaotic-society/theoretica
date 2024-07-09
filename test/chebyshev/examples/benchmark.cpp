
///
/// @file errors.h Example program for benchmarks.
///

#include "chebyshev.h"
#include <cmath>
using namespace ch;


double f(double x) {
	return x * std::sqrt(x);
}

double g(double x) {
	return std::atan(x * x);
}



int main(int argc, char const *argv[]) {

	// Setup benchmarking
	benchmark::setup("chebyshev", argc, argv);

		// Set options for multiple benchmarks
		benchmark::benchmark_options<double> opt {};
		opt.inputGenerator = benchmark::generator::uniform1D(0, 1000);
		opt.iterations = 1E+06;
		opt.runs = 10;

		// Benchmark the given functions
		benchmark::benchmark("f(x)", f, opt);
		benchmark::benchmark("g(x)", g, opt);

	// Stop benchmarking and exit
	benchmark::terminate();
}
