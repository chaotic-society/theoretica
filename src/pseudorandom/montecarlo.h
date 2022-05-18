
///
/// @file montecarlo.h Monte Carlo methods
///

#ifndef UROBORO_MONTECARLO_H
#define UROBORO_MONTECARLO_H

#include "./prng.h"
#include "./quasirandom.h"


namespace uroboro {


	/// Approximate an integral by using Monte Carlo integration
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param f_max The function maximum in the [a, b] interval
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	inline real approx_integral_montecarlo(
		real_function f,
		real a, real b, real f_max,
		PRNG& g, unsigned int N = 10000) {

		unsigned int N_inside = 0;

		for (int i = 0; i < N; ++i) {
		
			real x_n = rand_real(a, b, g);
			real y_n = rand_real(0, f_max, g);

			if(f(x_n) > y_n)
				N_inside++;
		}

		return (N_inside / static_cast<real>(N)) * (b - a) * f_max;
	}

}


#endif
