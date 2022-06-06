
///
/// @file benchmark_real_analysis.cpp Benchmark real functions
///

#ifndef THEORETICA_BENCHMARK_REAL_ANALYSIS
#define THEORETICA_BENCHMARK_REAL_ANALYSIS

#include "./timer.h"
#include "./benchmark.h"


int main(int argc, char const *argv[]) {

	// Setup
	setup_benchmark("real_analysis");
	print_benchmark_header();


	// Initialize random input data

	std::vector<real> input;
	input.resize(N);
	std::vector<real> input_norm;
	input_norm.resize(N);

	PRNG g = PRNG::xoshiro(
			std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::system_clock::now().time_since_epoch()).count());

	// Initialize input pool with random numbers in the interval [0, 1000000]
	// and input_norm in the interval [0, 1]
	for (unsigned int i = 0; i < N; ++i) {
		input[i] = rand_real(0, 1000000, g);
		input_norm[i] = rand_real(0, 1, g);
	}


	// Benchmark real functions

	benchmark_real_function("th::square", th::square, input);
	benchmark_real_function("th::cube", th::cube, input);
	benchmark_real_function("th::sqrt", th::sqrt, input);
	benchmark_real_function("th::cbrt", th::cbrt, input);
	benchmark_real_function("th::abs", th::abs, input);
	benchmark_real_function("th::fract", th::fract, input);
	benchmark_real_function("th::ln", th::ln, input);
	benchmark_real_function("th::log2", th::log2, input);
	benchmark_real_function("th::log10", th::log10, input);
	benchmark_real_function("th::atan", th::atan, input);
	benchmark_real_function("th::asin", th::asin, input_norm);
	benchmark_real_function("th::acos", th::acos, input_norm);
	benchmark_real_function("th::atan2", th::atan2, input, input);
	benchmark_real_function("th::exp", th::exp, input_norm);
	benchmark_real_function("th::powf", th::powf, input, input_norm);

	terminate_benchmark();
	
	return 0;
}


#endif
