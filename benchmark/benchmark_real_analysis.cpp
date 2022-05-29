
///
/// @file benchmark_real_analysis.cpp Benchmark real functions
///

#ifndef THEORETICA_BENCHMARK_REAL_ANALYSIS
#define THEORETICA_BENCHMARK_REAL_ANALYSIS

#include "./timer.h"
#include <iostream>

#include "theoretica.h"
using namespace th;


// Number of iterations
constexpr unsigned int N = 1000000;

// Number of runs
constexpr unsigned int M = 10;

// Name of the module being benchmarked
const char* module_name = "real_analysis";

// Current function
std::string curr_func_name = "";


void init_benchmark(const std::string& func_name) {
	curr_func_name = func_name;
	std::cout << "\nBenchmarking " << func_name << "..." << std::endl;
}


void end_benchmark(long double elapsed) {
	std::cout << "\n\tElapsed time: " << elapsed << " ms";
	std::cout << "\n\tComputations per second: "
			  << 1.0 / (elapsed * 0.001) << " /s\n" << std::endl;
	std::cout << "Ended benchmarking on " << curr_func_name << std::endl;
}


void benchmark_real_function(std::string func_name, real_function f,
							 const std::vector<real>& input) {

	init_benchmark(func_name);

	// Dummy variable
	__volatile__ real c = 0;

	// Sum of M runs with N iterations each
	long double sum = 0;

	for (int i = 0; i < M; ++i) {

		timer t = timer();

		for (int j = 0; j < N; ++j) {
			c += f(input[j]);
		}

		long double elapsed = t();
		sum += elapsed / (long double) N;
	}

	end_benchmark(sum / (long double) M);
}


int main(int argc, char const *argv[]) {

	std::cout << "Starting benchmark of real analysis" << std::endl;

	std::vector<real> input;
	input.resize(N);

	PRNG g = PRNG::xoshiro(
			std::chrono::duration_cast<std::chrono::seconds>(
			std::chrono::system_clock::now().time_since_epoch()).count());

	// Initialize input pool with random numbers in the interval [0, 1000000]
	for (int i = 0; i < N; ++i) {
		input[i] = rand_real(0, 1000000, g);
	}


	benchmark_real_function("th::square", th::sqrt, input);
	benchmark_real_function("th::cube", th::cube, input);
	benchmark_real_function("th::sqrt", th::sqrt, input);
	benchmark_real_function("th::cbrt", th::cbrt, input);
	benchmark_real_function("th::abs", th::abs, input);
	benchmark_real_function("th::fract", th::fract, input);
	benchmark_real_function("th::ln", th::ln, input);
	benchmark_real_function("th::log2", th::log2, input);
	benchmark_real_function("th::log10", th::log10, input);
	benchmark_real_function("th::atan", th::atan, input);
	benchmark_real_function("th::asin", th::asin, input);
	benchmark_real_function("th::acos", th::acos, input);
	
	return 0;
}


#endif
