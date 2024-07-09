
///
/// @file fail.h Default fail functions.
///

#ifndef CHEBYSHEV_FAIL_H
#define CHEBYSHEV_FAIL_H

#include "./prec_structures.h"


namespace chebyshev {
namespace prec {


	namespace fail {


		/// Default fail function which marks the test as failed
		/// if the maximum error on the domain is bigger than the tolerance
		auto fail_on_max_err = [](estimate_result r) -> bool {
			return (r.maxErr > r.tolerance) || (r.maxErr != r.maxErr);
		};


		/// Marks the test as failed if the mean error on the domain
		/// is bigger than the tolerance or the error is NaN.
		auto fail_on_mean_err = [](estimate_result r) -> bool {
			return (r.meanErr > r.tolerance) || (r.meanErr != r.meanErr);
		};


		/// Marks the test as failed if the RMS error on the domain
		/// is bigger than the tolerance or the error is NaN.
		auto fail_on_rms_err = [](estimate_result r) -> bool {
			return (r.rmsErr > r.tolerance) || (r.rmsErr != r.rmsErr);
		};


		/// Marks the test as failed if the relative error on the
		/// domain is bigger than the tolerance or the error is NaN.
		auto fail_on_rel_err = [](estimate_result r) -> bool {
			return (r.relErr > r.tolerance) || (r.relErr != r.relErr);
		};

	}
}}


#endif
