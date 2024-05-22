
///
/// @file benchmark.h Function benchmarking
///

#pragma once

/// Default number of iterations
#ifndef CHEBYSHEV_ITER
#define CHEBYSHEV_ITER 1000000
#endif

/// Default number of runs for benchmarks
#ifndef CHEBYSHEV_RUNS
#define CHEBYSHEV_RUNS 10
#endif
	
#include "core/common.h"
#include "core/timer.h"
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>


// Benchmark a real function on uniformly distributed values in [a, b]
#define BENCHMARK(f, a, b) chebyshev::benchmark::request(#f, [](chebyshev::Real x){ return f(x); }, chebyshev::uniform_generator(a, b))


namespace chebyshev {

	/// @namespace chebyshev::benchmark Function benchmarking
	namespace benchmark {

		/// Benchmark run request, used to store
		/// information about requested benchmarks
		/// for later execution
		struct benchmark_request {
			std::string funcName {"unknown"};
			RealFunction func {nullptr};
			RealInputGenerator gen {nullptr};
			unsigned int iter {CHEBYSHEV_ITER};
			unsigned int runs {CHEBYSHEV_RUNS};
		};


		/// Benchmark result, used to store
		/// information about a benchmark execution
		struct benchmark_result {
			std::string funcName {"unknown"};
			unsigned int iter {CHEBYSHEV_ITER};
			unsigned int runs {CHEBYSHEV_RUNS};
			Real total_time {0};
			Real avg_time {0};
			Real runs_per_sec {0};

			benchmark_result() {}

			benchmark_result(
				const std::string& name, Real elapsed_time,
				unsigned int i, unsigned int r) {

				funcName = name;
				total_time = elapsed_time;
				avg_time = elapsed_time / (i * r);
				runs_per_sec = 1.0 / (avg_time * 0.001);
				iter = i;
				runs = r;
			}


			benchmark_result(Real elapsed_time, unsigned int i, unsigned int r) {

				total_time = elapsed_time;
				avg_time = elapsed_time / (i * r);
				runs_per_sec = 1.0 / (avg_time * 0.001);
				iter = i;
				runs = r;
			}

		};


		/// A custom benchmark function which takes as input the number of iterations
		/// and runs requested and returns a benchmark result
		using CustomBenchmarkFunction = std::function<benchmark_result(unsigned int, unsigned int)>;


		/// A custom benchmark request
		struct benchmark_custom_request {
			std::string funcName {"unknown"};
			CustomBenchmarkFunction f;
			unsigned int iter {CHEBYSHEV_ITER};
			unsigned int runs {CHEBYSHEV_RUNS};
		};


		/// @class benchmark_state Global state of benchmarking
		struct benchmark_state {
			
			/// List of requested benchmark runs.
			/// Benchmarks are run  when the benchmark::run() function is called
			std::vector<benchmark_request> requests;

			/// List of requested custom benchmarks.
			std::vector<benchmark_custom_request> customRequests;

			/// Print to standard output?
			bool quiet = false;

			/// Name of the module currently being benchmarked
			std::string moduleName;

			/// Default number of iterations
			unsigned int defaultIterations = CHEBYSHEV_ITER;

			/// Default number of runs
			unsigned int defaultRuns = CHEBYSHEV_RUNS;

			/// Output file for the current module
			std::ofstream outputFile;

			/// Relative or absolute path to output folder
			std::string outputFolder = "";

			/// Benchmark state.results
			std::vector<benchmark_result> results;

			/// Number of failed benchmarks
			unsigned int failedBenchmarks = 0;

			/// Target benchmarks marked for execution
			/// (all benchmarks will be executed if empty)
			std::map<std::string, bool> pickedBenchmarks;
			
		} state;


		/// Setup a module's benchmark
		inline void setup(const std::string& module = "unknown",
			int argc = 0, const char** argv = nullptr,
			unsigned int iter = 0,
			unsigned int runs = 0) {

			// Initialize pick list
			if(argc && argv) {
				for (int i = 1; i < argc; ++i) {
					state.pickedBenchmarks[argv[i]] = true;
				}
			}

			state.moduleName = module;

			if(iter)
				state.defaultIterations = iter;

			if(runs)
				state.defaultRuns = runs;

			srand(time(nullptr));

			if(state.outputFile.is_open())
				state.outputFile.close();

			state.outputFile.open(state.outputFolder + "benchmark_" + state.moduleName + ".csv");

			if(!state.outputFile.is_open()) {
				std::cout << "Can't open output file" << std::endl;
				exit(1);
			}

		}


		/// Register a function to be benchmarked
		inline void request(const std::string& funcName,
			RealFunction f, RealInputGenerator g,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			benchmark_request r;
			r.funcName = funcName;
			r.func = f;
			r.gen = g;
			r.iter = n;
			r.runs = m;

			state.requests.push_back(r);
		}


		// /// Register a function to be benchmarked
		// void request(const std::string& funcName,
		// 	Real(*f)(Real), RealInputGenerator g,
		// 	unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

		// 	request(funcName, [f](Real x) {return f(x);}, g, n, m);
		// }


		/// Request a custom benchmark
		inline void custom_request(
			const std::string& funcName, CustomBenchmarkFunction f,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			benchmark_custom_request r;
			r.funcName = funcName;
			r.f = f;
			r.iter = n;
			r.runs = m;

			state.customRequests.push_back(r);
		}



		/// Benchmark a function
		inline benchmark_result benchmark(const std::string& funcName, RealFunction f, RealInputGenerator g,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			// Dummy variable
			__volatile__ Real c = 0;

			std::vector<Real> input;
			input.reserve(n);

			for (unsigned int i = 0; i < n; ++i)
				input[i] = g(i);

			// Sum of m runs with n iterations each
			long double sum = 0;

			for (unsigned int i = 0; i < m; ++i) {

				timer t = timer();

				for (unsigned int j = 0; j < n; ++j)
					c += f(input[j]);

				sum += t();
			}

			return benchmark_result(funcName, sum, n, m);
		}


		/// Benchmark a function
		inline benchmark_result benchmark(
			const std::string& funcName, RealFunction f,
			const std::vector<Real> input,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			if(input.size() < n) {
				std::cout << "Wrong input size in benchmark, skipping request ..." << std::endl;
				state.failedBenchmarks++;
				return benchmark_result();
			}

			// Dummy variable
			__volatile__ Real c = 0;

			// Sum of m runs n iterations each
			long double sum = 0;

			for (unsigned int i = 0; i < m; ++i) {

				timer t = timer();

				for (unsigned int j = 0; j < n; ++j)
					c += f(input[j]);

				sum += t();
			}

			return benchmark_result(funcName, sum, n, m);
		}


		/// Benchmark a function
		inline benchmark_result benchmark(benchmark_request r) {
			return benchmark(r.funcName, r.func, r.gen, r.iter, r.runs);
		}


		/// Print a benchmark result
		inline void print_benchmark(const benchmark_result& br) {
			std::cout << std::left << std::setw(20) << br.funcName << " | ";
			std::cout << std::setw(12) << br.avg_time << " | ";
			std::cout << std::setw(10) << std::right << std::floor(br.runs_per_sec) << std::endl;
		}


		/// Run all registered benchmarks
		inline void run() {

			std::cout.precision(8);

			if(!state.quiet) {
				std::cout << "Starting benchmark of " << state.moduleName << std::endl;
				std::cout << "Parameters: Iterations = " << state.defaultIterations << ", Runs = " << state.defaultRuns << "\n" << std::endl;

				std::cout << std::left << std::setw(20) << "Function" << " | "
				<< std::setw(12) << "Time (ms)" << " | " << std::setw(12) << "Runs/sec" << std::endl;
				state.outputFile << "Function, Time(ms), Runs/sec" << std::endl;
			}
		
			for (const auto& r : state.requests) {

				// Skip benchmark if it hasn't been picked
				if(!state.pickedBenchmarks.empty() && !state.pickedBenchmarks[r.funcName])
					continue;

				benchmark_result br = benchmark(r);
				state.results.push_back(br);
				
				if(!state.quiet)
					print_benchmark(br);

				state.outputFile << br.funcName << ", "
							<< br.avg_time << ", "
							<< br.runs_per_sec << std::endl;
			}

			for (const auto& r : state.customRequests) {

				// Skip benchmark if it hasn't been picked
				if(!state.pickedBenchmarks.empty() && !state.pickedBenchmarks[r.funcName])
					continue;

				benchmark_result br = r.f(r.iter, r.runs);
				br.funcName = r.funcName;

				state.results.push_back(br);
				
				if(!state.quiet)
					print_benchmark(br);

				state.outputFile << br.funcName << ", "
							<< br.avg_time << ", "
							<< br.runs_per_sec << std::endl;
			}

			state.requests.clear();
		}


		/// End benchmarking of the current module
		inline void terminate(bool exit = true) {

			if(state.requests.size() + state.customRequests.size())
				run();

			std::cout << "\nFinished benchmark of " << state.moduleName << std::endl;
			if(state.failedBenchmarks)
				std::cout << state.failedBenchmarks << " benchmarks failed!" << std::endl;

			std::string filename = state.outputFolder + "benchmark_" + state.moduleName + ".csv";
			std::cout << "Results have been saved in " << filename << std::endl;

			state.moduleName = "unknown";
			state.outputFile.close();

			if(exit)
				std::exit(state.failedBenchmarks);
		}

	}

}
