
///
/// @file prec_structures.h Structures for precision testing.
///

#ifndef CHEBYSHEV_PREC_STRUCTURES_H
#define CHEBYSHEV_PREC_STRUCTURES_H

#include <string>
#include <vector>
#include <functional>

#include "../core/common.h"
#include "./interval.h"
#include "./distance.h"


namespace chebyshev {

	namespace prec {


		/// @class estimate_result
		/// A structure holding the result of precision estimation.
		struct estimate_result {
			
			/// Identifying name of the function or test case.
			std::string name = "unknown";

			/// Interval of estimation.
			std::vector<interval> domain {};

			/// Tolerance on the max absolute error.
			long double tolerance = 0;

			/// Estimated maximum absolute error on interval.
			long double maxErr = get_nan<long double>();

			/// Estimated mean error on interval.
			long double meanErr = get_nan<long double>();

			/// Estimated RMS error on interval.
			long double rmsErr = get_nan<long double>();

			/// Estimated relative error on interval.
			long double relErr = get_nan<long double>();

			/// Estimated absolute error on interval.
			long double absErr = get_nan<long double>();

			/// Additional fields by name, as a floating point value.
			std::map<std::string, long double> additionalFields {};

			/// Whether the test failed.
			bool failed = false;

			/// Print to standard output or not.
			bool quiet = false;

			/// Total number of iterations for integral quadrature.
			unsigned int iterations {0};
		};


		/// A function which determines whether an estimation failed.
		using FailFunction = std::function<bool(const estimate_result&)>;


		/// Distance function between two elements.
		template<typename Type>
		using DistanceFunction = std::function<long double(Type, Type)>;


		/// @class estimate_options
		/// A structure holding the options for precision estimation.
		template<typename R, typename ...Args>
		struct estimate_options {

			using Estimator_t = std::function<
				estimate_result(
					std::function<R(Args...)>,
					std::function<R(Args...)>,
					estimate_options)>;

			/// The domain of estimation.
			std::vector<interval> domain {};

			/// The precision estimator to use
			/// (defaults to a dummy estimator)
			Estimator_t estimator = [](
				std::function<R(Args...)> f1,
				std::function<R(Args...)> f2,
				estimate_options opt) {
				return estimate_result();
			};

			/// The tolerance to use to determine whether the test failed.
			long double tolerance = CHEBYSHEV_PREC_TOLERANCE;

			/// Number of function evaluations to use.
			unsigned int iterations = CHEBYSHEV_PREC_ITER;

			/// The function to determine whether the test failed
			/// (defaults to fail::fail_on_max_err).
			FailFunction fail = [](const estimate_result& r) {
				return (r.maxErr > r.tolerance) || (r.maxErr != r.maxErr);
			};

			/// Whether to show the test result or not.
			bool quiet = false;


			/// Construct estimate options with all default values.
			/// @note The estimator and domain must be set to 
			/// correctly use the options for test cases.
			estimate_options() {}


			/// Construct estimate options from a one-dimensional
			/// interval domain and an estimator, with other fields
			/// equal to the default values.
			estimate_options(interval omega, Estimator_t estimator)
			: domain({omega}), estimator(estimator) {}


			/// Construct estimate options from a one-dimensional
			/// interval domain, an estimator, a tolerance and an optional
			/// quiet flag, with other fields equal to the default values.
			estimate_options(
				interval omega,
				Estimator_t estimator,
				long double tolerance,
				bool quiet = false)
			: domain({omega}), estimator(estimator), tolerance(tolerance), quiet(quiet) {}


			/// Construct estimate options from a multidimensional
			/// interval domain and an estimator, with other fields
			/// equal to the default values.
			estimate_options(std::vector<interval> omega, Estimator_t estimator)
			: domain(omega), estimator(estimator) {}

			/// Construct estimate options from a multidimensional
			/// interval domain, an estimator, a tolerance and an optional quiet flag,
			/// with other fields equal to the default values.
			estimate_options(
				std::vector<interval> omega,
				Estimator_t estimator,
				long double tolerance,
				bool quiet = false)
			: domain(omega), estimator(estimator), tolerance(tolerance), quiet(quiet) {}
			
		};


		/// Generic precision estimator function signature.
		template<typename R, typename ...Args>
		using Estimator = typename estimate_options<R, Args...>::Estimator_t;


		/// @class equation_result
		/// A structure holding the result of an evaluation.
		struct equation_result {

			/// Identifying name of the function or test case.
			std::string name = "unknown";

			/// Evaluated value.
			long double evaluated = get_nan<long double>();

			/// Expected value.
			long double expected = get_nan<long double>();
			
			/// Evaluated difference between expected and evaluated values.
			long double difference = get_nan<long double>();

			/// Additional fields by name,
			/// as a floating point value.
			std::map<std::string, long double> additionalFields {};

			/// Tolerance on the absolute difference.
			long double tolerance = 0;

			/// Whether the test failed.
			bool failed = true;

			/// Print to standard output or not.
			bool quiet = false;
		};


		/// @class equation_options
		/// Structure holding options for equivalence evaluation.
		template<typename T>
		struct equation_options {
			
			/// Tolerance on the absolute difference.
			long double tolerance = CHEBYSHEV_PREC_TOLERANCE;

			/// Distance function to measure the distance
			/// between the expected and evaluated value.
			DistanceFunction<T> distance = [](T x, T y) {
				const auto diff = x - y;
				return (long double) (diff > 0 ? diff : -diff);
			};

			/// Print to standard output or not.
			bool quiet = false;


			/// Default constructor for equation options.
			equation_options() {}


			/// Construct equation options from the tolerance,
			/// setting the distance function to a simple Euclidean distance.
			equation_options(long double tolerance) : tolerance(tolerance) {}


			/// Construct equation options from the tolerance,
			/// the distance function and the quiet flag (defaults to false).
			equation_options(long double tolerance, DistanceFunction<T> dist, bool quiet = false)
			: tolerance(tolerance), distance(dist), quiet(quiet) {}
		};

	}
}

#endif
