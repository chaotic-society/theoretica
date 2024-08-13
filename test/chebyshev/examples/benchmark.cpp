
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

unsigned int h(unsigned int n) {
	if(n == 0)
		return 0;
	else
		return n + h(n - 1);
}



int main(int argc, char const *argv[]) {

	// Setup benchmarking
	benchmark::setup("example", argc, argv);

		// Set the output file for the benchmark module
		benchmark::state.outputFiles = { "example_benchmark.csv" };

		// Set options for multiple benchmarks
		// with a benchmark_options structure,
		// specialized for functions taking in doubles
		auto opt = benchmark::benchmark_options<double>(
			10, 	// runs
			1E+06,	// iterations
			benchmark::generator::uniform1D(0, 1000) // input generator
		);

		// Benchmark the given functions
		benchmark::benchmark("f(x)", f, opt);
		benchmark::benchmark("g(x)", g, opt);

		// Specify parameters directly
		benchmark::benchmark<unsigned int>(
			"h(n)",
			h, 10, 1000,
			[](unsigned int i) {
				return random::natural() % 1000;
			}
		);

		// You may need to specify the input type
		// of your function if it isn't deduced.

	// Stop benchmarking and exit
	benchmark::terminate();
}
