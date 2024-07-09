
///
/// @file common.h Common definitions for the framework.
///

#ifndef CHEBYSHEV_COMMON_H
#define CHEBYSHEV_COMMON_H

#ifndef CHEBYSHEV_PREC_ITER
/// Default number of function evaluations
/// in precision testing.
#define CHEBYSHEV_PREC_ITER 1E+06
#endif

#ifndef CHEBYSHEV_PREC_TOLERANCE
/// Default tolerance in precision testing.
#define CHEBYSHEV_PREC_TOLERANCE 1E-08
#endif

#ifndef CHEBYSHEV_BENCHMARK_ITER
/// Default number of benchmark iterations.
#define CHEBYSHEV_BENCHMARK_ITER 1000
#endif

#ifndef CHEBYSHEV_BENCHMARK_RUNS
/// Default number of benchmark runs.
#define CHEBYSHEV_BENCHMARK_RUNS 10
#endif

#ifndef CHEBYSHEV_OUTPUT_WIDTH
/// Default width of output columns
#define CHEBYSHEV_OUTPUT_WIDTH 12
#endif


#include <limits>
#include <vector>
#include "../prec/interval.h"


namespace chebyshev {


	/// A real function of real variable.
	template<typename FloatType = double>
	using RealFunction = std::function<FloatType(FloatType)>;


	/// Get a quiet NaN of the specified floating point type.
	template<typename FloatType = long double>
	inline constexpr FloatType get_nan() {
		return std::numeric_limits<FloatType>::quiet_NaN();
	}


	/// Generate a pseudorandom number following
	/// a uniform distribution over the interval.
	inline long double random_uniform(long double a, long double b) {
		return (rand() / (long double) RAND_MAX) * (b - a) + a;
	}

	/// Generate a pseudorandom vector
	/// with uniformly distributed elements,
	/// overwriting the argument.
	template<typename Vector>
	inline void sample_uniform(Vector& x, const std::vector<prec::interval> domain) {

		if(x.size() != domain.size())
			throw std::runtime_error(
				"Vector and domain size mismatch in chebyshev::sample_uniform");

		for (int i = 0; i < x.size(); ++i)
			x[i] = random_uniform(domain[i].a, domain[i].b);
	}
	
}

#endif
