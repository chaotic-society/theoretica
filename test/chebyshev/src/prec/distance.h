
///
/// @file distance.h Distance functions.
///

#ifndef CHEBYSHEV_DISTANCE_H
#define CHEBYSHEV_DISTANCE_H


namespace chebyshev {
namespace prec {

	/// @namespace chebyshev::prec::distance Distance functions for use in prec::equals
	namespace distance {


		/// Absolute distance between two real values.
		template<typename FloatType = double>
		inline FloatType abs_distance(FloatType a, FloatType b) {

			const FloatType diff = b - a;
			return (diff > 0) ? diff : -diff;
		}


	}

}}

#endif
