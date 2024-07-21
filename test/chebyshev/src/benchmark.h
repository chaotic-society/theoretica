
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


	/// @namespace chebyshev::benchmark Benchmark module.
	namespace benchmark {


		/// @class benchmark_state
		/// Global state of the benchmark module.
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
			std::vector<std::string> benchmarkFiles {};

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
		///
		/// @param exit Whether to exit after terminating the module.
		inline void terminate(bool exit = true) {

			output::state.quiet = state.quiet;

			// Output to file is true but no specific files are specified,
			// add default output file.
			if(state.outputToFile && !state.benchmarkFiles.size()) {
				std::string filename;
				filename = output::state.outputFolder + state.moduleName + "_results";
				output::state.outputFiles[filename] = std::ofstream(filename);
			}

			output::print_results(state.benchmarkResults, state.benchmarkColumns, state.benchmarkFiles);

			std::cout << "Finished benchmarking " << state.moduleName << '\n';
			std::cout << state.totalBenchmarks << " total benchmarks, "
				<< state.failedBenchmarks << " failed (" << std::setprecision(3) << 
				(state.failedBenchmarks / (double) state.totalBenchmarks) * 100 << "%)"
				<< '\n';

			state = benchmark_state();

			if(exit) {
				output::terminate();
				std::exit(state.failedBenchmarks);
			}
		}


		/// Run a benchmark on a generic function,
		/// with the given options.
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& funcName,
			Function func,
			benchmark_options<InputType> opt) {

			// Generate input set
			std::vector<InputType> input (opt.iterations);
			for (unsigned int i = 0; i < opt.iterations; ++i)
				input[i] = opt.inputGenerator(i);

			// Sum of m runs with n iterations each
			long double totalRuntime = 0;

			// Whether the benchmark failed because of an exception
			bool failed = false;

			try {

				// Dummy variable
				__volatile__ auto c = func(input[0]);

				for (unsigned int i = 0; i < opt.runs; ++i) {

					timer t = timer();

					for (unsigned int j = 0; j < opt.iterations; ++j)
						c += func(input[j]);

					totalRuntime += t();
				}

				// Differentiate between types with operator+= or not ?
				// c = *((&c + 1) - 1);

			} catch(...) {
				failed = true;
			}

			benchmark_result res {};
			res.funcName = funcName;
			res.runs = opt.runs;
			res.iterations = opt.iterations;
			res.totalRuntime = totalRuntime;
			res.averageRuntime = totalRuntime / (opt.runs * opt.iterations);
			res.runsPerSecond = 1000.0 / res.averageRuntime;
			res.failed = failed;

			state.totalBenchmarks++;
			if(failed)
				state.failedBenchmarks++;

			state.benchmarkResults[funcName].push_back(res);
		}


		/// Run a benchmark on a generic function,
		/// with the given options.
		template<typename InputType = double, typename Function>
		inline void benchmark(
			const std::string& funcName,
			Function func,
			unsigned int runs = CHEBYSHEV_BENCHMARK_RUNS,
			unsigned int iterations = CHEBYSHEV_BENCHMARK_ITER,
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
