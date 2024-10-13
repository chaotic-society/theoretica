
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


		/// @class benchmark_settings Global settings of the benchmark module
		struct benchmark_settings {

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

			/// Target benchmarks marked for execution
			/// (all benchmarks will be executed if empty)
			std::map<std::string, bool> pickedBenchmarks {};

			/// The files to write benchmark results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> benchmarkOutputFiles {};

			/// Default columns to print for benchmarks.
			std::vector<std::string> benchmarkColumns = {
				"name", "averageRuntime", "stdevRuntime", "runsPerSecond"
			};
			
		} settings;


		/// @class benchmark_results Results of benchmarks.
		struct benchmark_results {
			
			/// Total number of benchmarks
			unsigned int totalBenchmarks = 0;

			/// Number of failed benchmarks
			unsigned int failedBenchmarks = 0;

			/// Results of the benchmarks.
			std::map<std::string, std::vector<benchmark_result>> benchmarkResults {};

		} results;


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
					settings.pickedBenchmarks[argv[i]] = true;

			std::cout << "Starting benchmarks of the "
				<< moduleName << " module ..." << std::endl;

			settings.moduleName = moduleName;
			results.totalBenchmarks = 0;
			results.failedBenchmarks = 0;

			random::setup();
			output::setup();
		}


		/// Terminate the benchmarking environment.
		/// If benchmarks have been run, their results will be printed.
		///
		/// @param exit Whether to exit after terminating the module.
		inline void terminate(bool exit = true) {

			output::settings.quiet = settings.quiet;

			// Output to file is true but no specific files are specified, add default output file.
			if(	 settings.outputToFile &&
				!output::settings.outputFiles.size() &&
				!settings.benchmarkOutputFiles.size() &&
				!settings.outputFiles.size()) {
				
				settings.outputFiles = { settings.moduleName + "_results" };
			}

			std::vector<std::string> outputFiles;

			// Print benchmark results
			outputFiles  = settings.outputFiles;
			outputFiles.insert(outputFiles.end(), settings.benchmarkOutputFiles.begin(), settings.benchmarkOutputFiles.end());

			output::print_results(results.benchmarkResults, settings.benchmarkColumns, outputFiles);

			std::cout << "Finished benchmarking " << settings.moduleName << '\n';
			std::cout << results.totalBenchmarks << " total benchmarks, "
				<< results.failedBenchmarks << " failed (" << std::setprecision(3) << 
				(results.failedBenchmarks / (double) results.totalBenchmarks) * 100 << "%)"
				<< '\n';

			// Discard previous results
			results = benchmark_results();

			if(exit) {
				output::terminate();
				std::exit(results.failedBenchmarks);
			}
		}


		/// Measure the total runtime of a function over
		/// the given input for a single run. It is generally
		/// not needed to call this function directly,
		/// as benchmarks can be run and registered using
		/// benchmark::benchmark.
		///
		/// @param func The function to measure the runtime of
		/// @param input The vector of inputs
		/// @return The total runtime of the function over the input vector.
		template<typename InputType, typename Function>
		inline long double runtime(Function func, const std::vector<InputType>& input) {

			if (input.size() == 0)
				return 0.0;

			// Dummy variable
			__volatile__ auto c = func(input[0]);

			timer t = timer();

			for (unsigned int j = 0; j < input.size(); ++j)
				c += func(input[j]);

			return t();
		}


		/// Run a benchmark on a generic function, with the given input vector.
		/// The result is registered inside results.benchmarkResults.
		///
		/// @param name The name of the test case
		/// @param func The function to benchmark
		/// @param input The vector of input values
		/// @param runs The number of runs with the same input
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& name,
			Function func,
			const std::vector<InputType>& input,
			unsigned int runs = settings.defaultRuns,
			bool quiet = false) {

			// Whether the benchmark failed because of an exception
			bool failed = false;

			// Running average
			long double averageRuntime;

			// Running total sum of squares
			long double sumSquares;

			// Total runtime
			long double totalRuntime;

			try {

				// Use Welford's algorithm to compute the average and the variance
				totalRuntime = runtime(func, input);
				averageRuntime = totalRuntime / input.size();
				sumSquares = 0.0;

				for (unsigned int i = 1; i < runs; ++i) {
					
					// Compute the runtime for a single run
					// and update the running estimates
					const long double currentRun = runtime(func, input);
					const long double currentAverage = currentRun / input.size();
					totalRuntime += currentRun;

					const long double tmp = averageRuntime;
					averageRuntime = tmp + (currentAverage - tmp) / (i + 1);
					sumSquares += (currentAverage - tmp)
						* (currentAverage - averageRuntime);
				}

			} catch(...) {

				// Catch any exception and mark the benchmark as failed
				failed = true;
			}

			benchmark_result res {};
			res.name = name;
			res.runs = runs;
			res.iterations = input.size();
			res.totalRuntime = totalRuntime;
			res.averageRuntime = averageRuntime;
			res.runsPerSecond = 1000.0 / res.averageRuntime;
			res.failed = failed;
			res.quiet = quiet;

			if (runs > 1)
				res.stdevRuntime = std::sqrt(sumSquares / (runs - 1));

			results.totalBenchmarks++;
			if(failed)
				results.failedBenchmarks++;

			results.benchmarkResults[name].push_back(res);
		}


		/// Run a benchmark on a generic function, with the given options.
		/// The result is registered inside results.benchmarkResults.
		///
		/// @param name The name of the test case
		/// @param func The function to benchmark
		/// @param opt The benchmark options
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& name,
			Function func,
			const benchmark_options<InputType>& opt) {

			// Generate input set
			std::vector<InputType> input (opt.iterations);
			for (unsigned int i = 0; i < opt.iterations; ++i)
				input[i] = opt.inputGenerator(i);

			// Benchmark over input set
			benchmark(name, func, input, opt.runs, opt.quiet);
		}


		/// Run a benchmark on a generic function, with the given argument options.
		/// The result is registered inside results.benchmarkResults.
		///
		/// @param name The name of the test case
		/// @param func The function to benchmark
		/// @param run The number of runs with the same input
		/// @param iterations The number of iterations of the function
		/// @param inputGenerator The input generator to use
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& name,
			Function func,
			unsigned int runs = settings.defaultRuns,
			unsigned int iterations = settings.defaultIterations,
			InputGenerator<InputType> inputGenerator = generator::uniform1D(0, 1),
			bool quiet = false) {

			benchmark_options<InputType> opt;
			opt.runs = runs;
			opt.iterations = iterations;
			opt.inputGenerator = inputGenerator;
			opt.quiet = quiet;

			benchmark(name, func, opt);
		}
	}
}

#endif
