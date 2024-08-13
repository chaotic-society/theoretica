
///
/// @file benchmark.h Benchmark module.
///

#ifndef CHEBYSHEV_BENCHMARK_H
#define CHEBYSHEV_BENCHMARK_H

#include <ctime>
#include <iostream>

#include "./core/random.h"
#include "./benchmark/timer.h"
#include "./benchmark/generator.h"
#include "./benchmark/benchmark_structures.h"


namespace chebyshev {


	/// @namespace chebyshev::benchmark Benchmark module
	///
	/// This module provides routines for measuring the average
	/// runtime of functions of any kind over a randomized or fixed
	/// vector of inputs. The benchmark::benchmark implements this
	/// functionality and registers the results for analysis and output.
	namespace benchmark {


		/// @class benchmark_state Global state of the benchmark module
		struct benchmark_state {

			/// Whether to print benchmark results
			/// to standard output.
			bool quiet = false;

			/// Name of the module currently being benchmarked
			std::string moduleName;

			/// Default number of iterations
			unsigned int defaultIterations = CHEBYSHEV_BENCHMARK_ITER;

			/// Default number of runs
			unsigned int defaultRuns = CHEBYSHEV_BENCHMARK_RUNS;

			/// Whether to output results to a file.
			bool outputToFile = true;

			/// The files to write all benchmark results to.
			std::vector<std::string> outputFiles {};

			/// Benchmark state.results
			std::vector<benchmark_result> results;

			/// Total number of benchmarks
			unsigned int totalBenchmarks = 0;

			/// Number of failed benchmarks
			unsigned int failedBenchmarks = 0;

			/// Target benchmarks marked for execution
			/// (all benchmarks will be executed if empty)
			std::map<std::string, bool> pickedBenchmarks {};

			/// The files to write benchmark results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> benchmarkOutputFiles {};

			/// Results of the benchmarks.
			std::map<std::string, std::vector<benchmark_result>> benchmarkResults {};

			/// Default columns to print for benchmarks.
			std::vector<std::string> benchmarkColumns = {
				"funcName", "averageRuntime", "runsPerSecond"
			};
			
		} state;


		/// Setup the benchmark environment.
		///
		/// @param moduleName Name of the module under test.
		/// @param argc The number of command line arguments
		/// @param argv A list of C-style strings containing
		/// the command line arguments.
		inline void setup(
			std::string moduleName,
			int argc = 0,
			const char** argv = nullptr) {

			// Initialize list of picked tests
			if(argc && argv)
				for (int i = 1; i < argc; ++i)
					state.pickedBenchmarks[argv[i]] = true;

			std::cout << "Starting benchmarks of the "
				<< moduleName << " module ..." << std::endl;

			state.moduleName = moduleName;
			state.failedBenchmarks = 0;

			random::setup();
			output::setup();
		}


		/// Terminate the benchmarking environment.
		/// If benchmarks have been run, their results will be printed.
		///
		/// @param exit Whether to exit after terminating the module.
		inline void terminate(bool exit = true) {

			output::state.quiet = state.quiet;

			// Output to file is true but no specific files are specified, add default output file.
			if(	 state.outputToFile &&
				!state.benchmarkOutputFiles.size() &&
				!state.outputFiles.size()) {
				
				state.outputFiles = { state.moduleName + "_results" };
			}

			std::vector<std::string> outputFiles;

			// Print benchmark results
			outputFiles  = state.outputFiles;
			outputFiles.insert(outputFiles.end(), state.benchmarkOutputFiles.begin(), state.benchmarkOutputFiles.end());

			output::print_results(state.benchmarkResults, state.benchmarkColumns, outputFiles);

			std::cout << "Finished benchmarking " << state.moduleName << '\n';
			std::cout << state.totalBenchmarks << " total benchmarks, "
				<< state.failedBenchmarks << " failed (" << std::setprecision(3) << 
				(state.failedBenchmarks / (double) state.totalBenchmarks) * 100 << "%)"
				<< '\n';

			// Reset module information
			state = benchmark_state();

			if(exit) {
				output::terminate();
				std::exit(state.failedBenchmarks);
			}
		}


		/// Measure the total runtime of a function over
		/// the given input for many runs. It is generally
		/// not needed to call this function directly,
		/// as benchmarks can be run and registered using
		/// benchmark::benchmark.
		///
		/// @param func The function to measure the runtime of
		/// @param input The vector of inputs
		/// @param runs The number of runs to make with the same input
		/// @return The total runtime of the function over the input
		/// vector and over many runs.
		template<typename InputType, typename Function>
		inline long double runtime(
			Function func,
			const std::vector<InputType>& input,
			unsigned int runs = state.defaultRuns) {

			if (input.size() == 0)
				return 0.0;

			long double totalRuntime = 0.0;

			// Dummy variable
			__volatile__ auto c = func(input[0]);

			for (unsigned int i = 0; i < runs; ++i) {

				timer t = timer();

				for (unsigned int j = 0; j < input.size(); ++j)
					c += func(input[j]);

				totalRuntime += t();
			}

			return totalRuntime;
		}


		/// Run a benchmark on a generic function, with the given input vector.
		/// The result is registered inside state.benchmarkResults.
		///
		/// @param funcName The name of the test case
		/// @param func The function to benchmark
		/// @param input The vector of input values
		/// @param runs The number of runs with the same input
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& funcName,
			Function func,
			const std::vector<InputType>& input,
			unsigned int runs = state.defaultRuns) {

			// Sum of m runs with n iterations each
			long double totalRuntime = get_nan();

			// Whether the benchmark failed because of an exception
			bool failed = false;

			try {

				// Measure the total runtime
				totalRuntime = runtime(func, input, runs);

			} catch(...) {

				// Catch any exception and mark the benchmark as failed
				failed = true;
			}

			benchmark_result res {};
			res.funcName = funcName;
			res.runs = runs;
			res.iterations = input.size();
			res.totalRuntime = totalRuntime;
			res.averageRuntime = totalRuntime / (runs * input.size());
			res.runsPerSecond = 1000.0 / res.averageRuntime;
			res.failed = failed;

			state.totalBenchmarks++;
			if(failed)
				state.failedBenchmarks++;

			state.benchmarkResults[funcName].push_back(res);
		}


		/// Run a benchmark on a generic function, with the given options.
		/// The result is registered inside state.benchmarkResults.
		///
		/// @param funcName The name of the test case
		/// @param func The function to benchmark
		/// @param opt The benchmark options
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& funcName,
			Function func,
			const benchmark_options<InputType>& opt) {

			// Generate input set
			std::vector<InputType> input (opt.iterations);
			for (unsigned int i = 0; i < opt.iterations; ++i)
				input[i] = opt.inputGenerator(i);

			// Benchmark over input set
			benchmark(funcName, func, input, opt.runs);
		}


		/// Run a benchmark on a generic function, with the given argument options.
		/// The result is registered inside state.benchmarkResults.
		///
		/// @param funcName The name of the test case
		/// @param func The function to benchmark
		/// @param run The number of runs with the same input
		/// @param iterations The number of iterations of the function
		/// @param inputGenerator The input generator to use
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& funcName,
			Function func,
			unsigned int runs = state.defaultRuns,
			unsigned int iterations = state.defaultIterations,
			InputGenerator<InputType> inputGenerator = generator::uniform1D(0, 1)) {

			benchmark_options<InputType> opt;
			opt.runs = runs;
			opt.iterations = iterations;
			opt.inputGenerator = inputGenerator;

			benchmark(funcName, func, opt);
		}
	}
}

#endif
