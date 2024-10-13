
///
/// @file benchmark_structures.h Structures for the benchmark module.
///

#ifndef CHEBYSHEV_BENCHMARK_STRUCTURES_H
#define CHEBYSHEV_BENCHMARK_STRUCTURES_H

#include <functional>
#include <map>

#include "../core/common.h"
#include "./generator.h"


namespace chebyshev {

	namespace benchmark {


		/// @class benchmark_result
		/// Structure holding the results of a benchmark.
		struct benchmark_result {
			
			/// Identifying name of the function or test case.
			std::string name = "unknown";

			/// Number of runs.
			unsigned int runs = 0;

			/// Number of iterations.
			unsigned int iterations = 0;

			/// Total runtime over all runs and iterations.
			long double totalRuntime = get_nan<long double>();

			/// Estimated average runtime.
			long double averageRuntime = get_nan<long double>();

			/// Sample standard deviation of the runtime.
			long double stdevRuntime = get_nan<long double>();

			/// Number of runs per second.
			long double runsPerSecond = get_nan<long double>();

			/// Whether the benchmark failed because
			/// an exception was thrown.
			bool failed = true;

			/// Whether to print to standard output or not.
			bool quiet = false;

			/// Additional fields in floating point representation.
			std::map<std::string, long double> additionalFields {};

		};


		/// A function which takes in an index and returns
		/// a generated input element.
		template<typename InputType>
		using InputGenerator = std::function<InputType(unsigned int)>;


		/// @class benchmark_options
		/// A structure holding the options of a benchmark.
		template<typename InputType = double>
		struct benchmark_options {
			
			/// Number of runs (run with the same input values).
			unsigned int runs = CHEBYSHEV_BENCHMARK_RUNS;

			/// Number of iterations.
			unsigned int iterations = CHEBYSHEV_BENCHMARK_ITER;

			/// The function to use to generate input for the benchmark.
			InputGenerator<InputType> inputGenerator = generator::uniform1D(0, 1);

			/// Whether to print to standard output or not.
			bool quiet = false;


			/// Default constructor for benchmark options.
			benchmark_options() {}

			/// Construct benchmark options from the number of runs and iterations
			/// and whether to print the case to output (defaults to false).
			benchmark_options(unsigned int runs, unsigned int iterations, bool quiet = false)
			: runs(runs), iterations(iterations), quiet(quiet) {}

			/// Construct benchmark options from the number of runs and iterations,
			/// the input generator to use and whether to write the result to output.
			benchmark_options(
				unsigned int runs,
				unsigned int iterations,
				InputGenerator<InputType> gen,
				bool quiet = false)
			: runs(runs), iterations(iterations), inputGenerator(gen), quiet(quiet) {}

		};

	}
}

#endif
