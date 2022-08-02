
///
/// @file definitions.h General definitions
///

#pragma once

#include "interval.h"
#include <functional>
#define _hypot(a, b) hypot(a, b)
#include <cmath>


/// Construct a RealFunction from any function
#define REAL_LAMBDA(f) [](chebyshev::Real x){ return f(x); }


namespace chebyshev {


	/// Real number type (defaults to double)
#ifdef CHEBYSHEV_FLOAT
	using Real = float;
#elif defined(CHEBYSHEV_LONG_DOUBLE)
	using Real = long double;
#else
	using Real = double;
#endif


	/// A real function of real argument
	using RealFunction = std::function<Real(Real)>;


	/// An input generating function
	using RealInputGenerator = std::function<Real(unsigned int)>;


	/// Returns a real random number generator which
	/// generates uniform numbers inside the interval k
	RealInputGenerator uniform_generator(interval k) {
		return [k](unsigned int i) {
			return (rand() / (Real) RAND_MAX) * k.length() + k.a;
		};
	}


	/// Returns a real random number generator which
	/// generates uniform numbers inside the interval [a, b]
	RealInputGenerator uniform_generator(Real a, Real b) {
		return [a, b](unsigned int i) {
			return (rand() / (Real) RAND_MAX) * (b - a) + a;
		};
	}

}
