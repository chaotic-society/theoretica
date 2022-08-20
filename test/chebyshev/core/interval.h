
///
/// @file interval.h An interval on the real numbers
///

#pragma once


namespace chebyshev {

	/// @class interval An interval on the real numbers
	struct interval {

		// [a, b]

		long double a;
		long double b;

		interval() : a(0), b(1) {}


		/// Construct an interval from its left and right bounds
		interval(long double l, long double r) : a(l), b(r) {}


		/// Returns the length of the interval
		inline long double length() const {
			const long double diff = b - a;
			return diff > 0 ? diff : -diff;
		}
	};

}
