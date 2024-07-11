
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
