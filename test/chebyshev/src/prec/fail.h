
///
/// @file fail.h Default fail functions.
///

#ifndef CHEBYSHEV_FAIL_H
#define CHEBYSHEV_FAIL_H

#include "./prec_structures.h"


namespace chebyshev {
namespace prec {

	/// @namespace chebyshev::prec::fail Fail functions for use in prec::estimate
	///
	/// Fail functions are used to evaluate whether a certain test
	/// case has failed error estimation. If the fail function returns
	/// true, the test has failed. For example, a test case may fail
	/// if the maximum error is over a certain threshold.
	namespace fail {


		/// Passthrough fail function which marks all tests as passed (not failed).
		inline auto passthrough() {
			return [](const estimate_result& r) -> bool {
				return false;
			};
		} 


		/// Default fail function which marks the test as failed
		/// if the maximum error on the domain is bigger than the tolerance
		inline auto fail_on_max_err() {
			return [](const estimate_result& r) -> bool {
				return (r.maxErr > r.tolerance) || (r.maxErr != r.maxErr);
			};
		}


		/// Marks the test as failed if the mean error on the domain
		/// is bigger than the tolerance or the error is NaN.
		inline auto fail_on_mean_err() {
			return [](const estimate_result& r) -> bool {
				return (r.meanErr > r.tolerance) || (r.meanErr != r.meanErr);
			};
		}


		/// Marks the test as failed if the RMS error on the domain
		/// is bigger than the tolerance or the error is NaN.
		inline auto fail_on_rms_err() {
			return [](const estimate_result& r) -> bool {
				return (r.rmsErr > r.tolerance) || (r.rmsErr != r.rmsErr);
			};
		}


		/// Marks the test as failed if the relative error on the
		/// domain is bigger than the tolerance or the error is NaN.
		inline auto fail_on_rel_err() {
			return [](const estimate_result& r) -> bool {
				return (r.relErr > r.tolerance) || (r.relErr != r.relErr);
			};
		}

	}
}}


#endif
