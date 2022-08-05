
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
		};


		/// @class benchmark_state Global state of benchmarking
		struct benchmark_state {
			
			/// List of requested benchmark runs.
			/// Benchmarks are run in random order
			/// when the benchmark::run() function is called
			std::vector<benchmark_request> requests;

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

		} state;


		/// Setup a module's benchmark
		void setup(const std::string& module = "unknown",
			unsigned int iter = CHEBYSHEV_ITER,
			unsigned int runs = CHEBYSHEV_RUNS) {

			state.moduleName = module;
			state.defaultIterations = iter;
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
		void request(const std::string& f_name,
			RealFunction f, RealInputGenerator g,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			benchmark_request r;
			r.funcName = f_name;
			r.func = f;
			r.gen = g;
			r.iter = n;
			r.runs = m;

			state.requests.push_back(r);
		}


		// /// Register a function to be benchmarked
		// void request(const std::string& f_name,
		// 	Real(*f)(Real), RealInputGenerator g,
		// 	unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

		// 	request(f_name, [f](Real x) {return f(x);}, g, n, m);
		// }


		/// Benchmark a function
		benchmark_result benchmark(const std::string& f_name, RealFunction f, RealInputGenerator g,
			unsigned int n = state.defaultIterations, unsigned int m = state.defaultRuns) {

			// Dummy variable
			__volatile__ Real c = 0;

			std::vector<Real> input;
			input.resize(n);

			for (unsigned int i = 0; i < n; ++i)
				input[i] = g(i);

			// Sum of m runs with n iterations each
			long double sum = 0;

			for (unsigned int i = 0; i < m; ++i) {

				timer t = timer();

				for (unsigned int j = 0; j < n; ++j)
					c += f(input[j]);

				long double elapsed = t();
				sum += elapsed / (long double) n;
			}

			benchmark_result br;
			br.funcName = f_name;
			br.iter = n;
			br.runs = m;
			br.total_time = sum * n;
			br.avg_time = sum / (long double) m;
			br.runs_per_sec = 1.0 / (sum * 0.001 / (long double) m);

			return br;
		}


		/// Benchmark a function
		benchmark_result benchmark(
			const std::string& f_name, RealFunction f,
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

				long double elapsed = t();
				sum += elapsed / (long double) n;
			}

			benchmark_result br;
			br.funcName = f_name;
			br.iter = n;
			br.runs = m;
			br.total_time = sum * n;
			br.avg_time = sum / (long double) m;
			br.runs_per_sec = 1.0 / (sum * 0.001 / (long double) m);

			return br;
		}


		/// Benchmark a function
		benchmark_result benchmark(benchmark_request r) {
			return benchmark(r.funcName, r.func, r.gen, r.iter, r.runs);
		}


		/// Run all registered benchmarks
		void run() {

			std::cout.precision(8);

			if(!state.quiet) {
				std::cout << "Starting benchmark of " << state.moduleName << std::endl;
				std::cout << "Parameters: Iterations = " << state.defaultIterations << ", Runs = " << state.defaultRuns << "\n" << std::endl;

				std::cout << std::left << std::setw(20) << "Function" << " | "
				<< std::setw(12) << "Time (ms)" << " | " << std::setw(12) << "Runs/sec" << std::endl;
				state.outputFile << "Function, Time(ms), Runs/sec" << std::endl;
			}
		
			for (const auto& r : state.requests) {

				benchmark_result br = benchmark(r);
				state.results.push_back(br);
				
				std::cout << std::left << std::setw(20) << br.funcName << " | "
				<< std::setw(12) << br.avg_time << " | "
				<< std::setw(10) << std::right << std::floor(br.runs_per_sec) << std::endl;

				state.outputFile << br.funcName << ", "
							<< br.avg_time << ", "
							<< br.runs_per_sec << std::endl;
			}

			state.requests.clear();
		}


		/// End benchmarking of the current module
		void terminate(bool exit = true) {

			if(state.requests.size())
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
