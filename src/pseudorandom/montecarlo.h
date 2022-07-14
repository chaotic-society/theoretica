
///
/// @file montecarlo.h Monte Carlo methods
///

#ifndef THEORETICA_MONTECARLO_H
#define THEORETICA_MONTECARLO_H

#include "./prng.h"
#include "./quasirandom.h"


namespace theoretica {


	/// Approximate an integral by using Crude Monte Carlo integration
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	inline real approx_integral_crude(
		real_function f,
		real a, real b,
		PRNG& g, unsigned int N = 10000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i)			
			sum_y += f(rand_uniform(a, b, g));

		return (b - a) * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Hit-or-miss Monte Carlo integration
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param f_max The function maximum in the [a, b] interval
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	inline real approx_integral_hom(
		real_function f,
		real a, real b, real f_max,
		PRNG& g, unsigned int N = 10000) {

		unsigned int N_inside = 0;

		for (unsigned int i = 0; i < N; ++i) {
		
			real x_n = rand_uniform(a, b, g);
			real y_n = rand_uniform(0, f_max, g);

			if(f(x_n) > y_n)
				N_inside++;
		}

		return (N_inside / static_cast<real>(N)) * (b - a) * f_max;
	}

}


#endif
