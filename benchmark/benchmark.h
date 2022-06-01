
///
/// @file benchmark.h Benchmarking functions
///

#ifndef THEORETICA_BENCHMARK_H
#define THEORETICA_BENCHMARK_H

#include "theoretica.h"
using namespace th;

#include <iostream>
#include <iomanip>
#include <fstream>

// Number of iterations
unsigned int N = 1000000;

// Number of runs
unsigned int M = 10;

// Name of the module being benchmarked
std::string module_name = "unknown";

// Current function
std::string curr_func_name = "unknown";

// Output file
std::ofstream output_file;


void setup_benchmark(std::string module, unsigned int iter = 1000000, unsigned int runs = 10) {

	module_name = module;
	N = iter;
	M = runs;
	output_file.open(std::string("./benchmark/") + module_name + ".csv");
}


void print_benchmark_header() {

	std::cout << "Starting benchmark of " << module_name << std::endl;
	std::cout << "Parameters: N = " << N << ", M = " << M << std::endl;
	
	for (int i = 0; i < 80; ++i)
		std::cout << '-';
	std::cout << std::endl;

	std::cout << std::left << "Function\t\tTime (ms)\tRuns/sec" << std::endl;
	output_file << "Function, Time(ms), Runs/sec" << std::endl;

	for (int i = 0; i < 80; ++i)
		std::cout << '-';
	std::cout << std::endl;
}


void terminate_benchmark() {
	std::cout << "\nFinished benchmark of " << module_name << std::endl;

	output_file.close();
}


// Initialize benchmark of a specific function
void init_benchmark(const std::string& func_name) {
	curr_func_name = func_name;
	std::cout << std::setw(12) << func_name << "\t\t";
}


// End benchmark of a specific function
void end_benchmark(long double elapsed) {

	unsigned long int comp_sec = floor(1.0 / (elapsed * 0.001));

	std::cout << std::setw(8) << elapsed
			  << "\t" << std::setw(8) << comp_sec << std::endl;

	output_file << curr_func_name << ", "
				<< elapsed << ", "
				<< comp_sec << std::endl;
	
	curr_func_name = "unknown";
}


// Automatically benchmark a real function
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


// Automatically benchmark a real function of two real parameters
void benchmark_real_function(std::string func_name, real(*f)(real, real),
							 const std::vector<real>& input1,
							 const std::vector<real>& input2) {

	init_benchmark(func_name);

	// Dummy variable
	__volatile__ real c = 0;

	// Sum of M runs with N iterations each
	long double sum = 0;

	for (int i = 0; i < M; ++i) {

		timer t = timer();

		for (int j = 0; j < N; ++j) {
			c += f(input1[j], input2[j]);
		}

		long double elapsed = t();
		sum += elapsed / (long double) N;
	}

	end_benchmark(sum / (long double) M);
}


#endif
