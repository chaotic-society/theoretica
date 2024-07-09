
///
/// @file benchmark.h Benchmark module.
///

#ifndef CHEBYSHEV_BENCHMARK_H
#define CHEBYSHEV_BENCHMARK_H

#include <ctime>

#include "benchmark/timer.h"
#include "benchmark/generator.h"
#include "benchmark/benchmark_structures.h"


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

			/// Output file for the current module
			std::ofstream outputFile;

			/// Relative or absolute path to output folder
			std::string outputFolder = "";

			/// Prefix to prepend to the filename, in addition
			/// to the module name.
			std::string filenamePrefix = "benchmark_";

			/// Suffix to append to the filename, in addition
			/// to the module name.
			std::string filenameSuffix = ".csv";

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

			/// Results of the benchmarks.
			std::map<std::string, std::vector<benchmark_result>> benchmarkResults {};
			
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

			srand(time(nullptr));
			output::setup();
		}


		/// Terminate the precision testing environment.
		///
		/// @param exit Whether to exit after terminating testing.
		inline void terminate(bool exit = true) {


			if(state.outputToFile) {

				std::string filename;
				filename = state.outputFolder + state.filenamePrefix
					+ state.moduleName + state.filenameSuffix;

				if(state.outputFile.is_open())
					state.outputFile.close();

				state.outputFile.open(filename);

				if(!state.outputFile.is_open()) {
					std::cout << "Unable to open output file,"
						" results will NOT be saved!" << std::endl;
					state.outputToFile = false;
				}
			}


			output::table_state benchmarkTable {};

			if(state.benchmarkResults.size()) {

				if(!state.quiet) {
					std::cout << "\n";
					output::header_benchmark(benchmarkTable);
				}

				// Print to file as CSV
				if(state.outputToFile)
					output::header_benchmark(benchmarkTable, state.outputFile);
			}


			for (auto it = state.benchmarkResults.begin();
				it != state.benchmarkResults.end(); ++it) {

				const auto res_list = it->second;

				for (size_t i = 0; i < res_list.size(); ++i) {

					benchmarkTable.rowIndex++;

					if(it != state.benchmarkResults.end()
					&& std::next(it) == state.benchmarkResults.end()
					&& (i == res_list.size() - 1))
						benchmarkTable.isLastRow = true;

					if(!state.quiet)
						output::print_benchmark(res_list[i], benchmarkTable);
				
					if(state.outputToFile)
						output::print_benchmark(
							res_list[i], benchmarkTable, state.outputFile);
				}
			}

			std::cout << "\nFinished benchmarking " << state.moduleName << '\n';
			std::cout << state.totalBenchmarks << " total benchmarks, "
				<< state.failedBenchmarks << " failed (" <<
				(state.failedBenchmarks / (double) state.totalBenchmarks) * 100 << "%)"
				<< '\n';
				
			std::cout << "Results have been saved in "
				<< state.outputFolder << state.filenamePrefix
				<< state.moduleName << state.filenameSuffix << std::endl;

			if(state.outputFile.is_open())
				state.outputFile.close();

			state = benchmark_state();

			if(exit)
				std::exit(state.failedBenchmarks);
		}


		/// Run a benchmark on a generic function,
		/// with the given options.
		template<typename Function, typename InputType = double>
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
		template<typename Function, typename InputType = double>
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
