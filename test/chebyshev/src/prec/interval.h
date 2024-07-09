
///
/// @file interval.h Interval over the real numbers.
///

#ifndef CHEBYSHEV_INTERVAL_H
#define CHEBYSHEV_INTERVAL_H


namespace chebyshev {

	namespace prec {

		/// @class interval
		/// An interval on the real numbers.
		struct interval {

			/// Lower extreme of the interval.
			long double a;

			/// Upper extreme of the interval.
			long double b;

			/// Construct over the origin.
			interval() : a(0), b(0) {}


			/// Construct an interval from its lower and upper bounds.
			interval(long double l, long double r) : a(l), b(r) {}


			/// Returns the length of the interval
			inline long double length() const {
				const long double diff = b - a;
				return diff > 0 ? diff : -diff;
			}
		};

	}
}


#endif
