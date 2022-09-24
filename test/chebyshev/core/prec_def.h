#pragma once

#include "common.h"
#include "interval.h"
#include <vector>


namespace chebyshev {
	
	namespace prec {

		/// @class estimate_result The result of error estimation
		struct estimate_result {
			
			/// Uniquely identifying name of the function
			std::string funcName = "unknown";

			/// Interval of estimation
			interval k;

			/// Tolerance on the max absolute error
			Real tolerance;

			/// Estimated maximum absolute error on interval
			Real max_err;

			/// Estimated mean error on interval
			Real mean_err;

			/// Estimated RMS error on interval
			Real rms_err;

			/// Estimated relative error on interval
			Real rel_err;

			/// Estimated absolute error on interval
			Real abs_err;

			/// Did the test fail?
			bool failed = false;

			/// Print to standard output or not
			bool quiet = false;

			/// Total number of iterations for integral quadrature
			uint32_t iterations;
		};


		/// @class equation_result The result of equation checking
		struct equation_result {

			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// Evaluated value
			Real evaluated;

			/// Expected value
			Real expected;
			
			/// Absolute difference between expected and evaluated values
			Real diff;

			/// Tolerance on the absolute difference
			Real tolerance;

			/// Did the test fail?
			bool failed;

			/// Print to standard output or not
			bool quiet = false;
		};


		/// A function which determines whether an estimation failed
		using FailFunction = std::function<bool(estimate_result)>;


		/// Default fail function which marks the test as failed
		/// if the maximum error on the domain is bigger than the tolerance
		FailFunction fail_on_max_err = [](estimate_result r) {
			return r.max_err > r.tolerance;
		};

		/// Marks the test as failed if the mean error on the domain is bigger than the tolerance
		FailFunction fail_on_mean_err = [](estimate_result r) {
			return r.mean_err > r.tolerance;
		};

		/// Marks the test as failed if the RMS error on the domain is bigger than the tolerance
		FailFunction fail_on_rms_err = [](estimate_result r) {
			return r.rms_err > r.tolerance;
		};


		/// Marks the test as failed if the relative error on the domain is bigger than the tolerance
		FailFunction fail_on_rel_err = [](estimate_result r) {
			return r.rel_err > r.tolerance;
		};


		/// @class estimate_request A precision estimation request
		struct estimate_request {
			
			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// The function to estimate
			RealFunction func = nullptr;

			/// A function returning the expected output
			RealFunction funcExpected = nullptr;

			/// Requested estimation intervals
			std::vector<interval> intervals;

			/// Precision testing tolerance on max absolute error
			Real tolerance;

			/// Number of iterations for integral quadrature
			uint32_t iterations;

			/// Fail function which determines whether the test failed
			FailFunction fail;

			/// Print to standard output or not
			bool quiet = false;
		};

		/// A custom precision estimation function taking as input a list of
		/// the requested intervals, the tolerance and the number of iterations
		using CustomEstimateFunction
		= std::function<estimate_result(interval, Real, uint32_t)>;


		/// @class estimate_request A precision estimation request
		struct estimate_custom_request {
			
			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// A custom precision estimation function
			CustomEstimateFunction f;

			/// Requested estimation intervals
			std::vector<interval> intervals;

			/// Precision testing tolerance on max absolute error
			Real tolerance;

			/// Number of iterations for integral quadrature
			uint32_t iterations;

			/// Print to standard output or not
			bool quiet = false;
		};


		/// @class equation_request An equation request
		struct equation_request {

			/// Uniquely identifying function name
			std::string funcName = "unknown";

			/// Evaluated value
			Real evaluated;

			/// Expected value
			Real expected;

			/// Tolerance
			Real tolerance;

			/// Print to standard output or not
			bool quiet = false;
		};

	}
}
