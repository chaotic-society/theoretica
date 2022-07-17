
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
	inline real integral_crude(
		real_function f,
		real a, real b,
		PRNG& g, unsigned int N = 10000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i)			
			sum_y += f(rand_uniform(a, b, g));

		return (b - a) * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Crude
	/// Quasi-Monte Carlo integration by sampling
	/// from the Weyl sequence
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param N The number of points to generate
	inline real integral_quasi_crude(
		real_function f,
		real a, real b,
		unsigned int N = 10000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i)			
			sum_y += f(a + qrand_weyl(i + 1) * (b - a));

		return (b - a) * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Hit-or-miss Monte Carlo integration
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param f_max The function maximum in the [a, b] interval
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	inline real integral_hom(
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


	/// Approximate an integral by using Hit-or-miss
	/// Quasi-Monte Carlo integration, sampling points
	/// from the Weyl bi-dimensional sequence
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param f_max The function maximum in the [a, b] interval
	/// @param N The number of points to generate
	inline real integral_quasi_hom(
		real_function f,
		real a, real b, real f_max,
		unsigned int N = 10000) {

		unsigned int N_inside = 0;

		for (unsigned int i = 0; i < N; ++i) {

			vec2 v = qrand_weyl2(i + 1);
		
			const real x_n = a + (b - a) * v[0];
			const real y_n = v[1] * f_max;

			if(f(x_n) > y_n)
				N_inside++;
		}

		return (N_inside / static_cast<real>(N)) * (b - a) * f_max;
	}


	/// Use the Hit-or-Miss Monte Carlo method
	/// to approximate a double integral
	/// @param f The multivariate function to integrate
	/// @param a The lower extreme of integration on x
	/// @param b The upper extreme of integration on x
	/// @param c The lower extreme of integration on y
	/// @param d The upper extreme of integration on y
	/// @param f_max The function maximum in the [a, b]x[c, d] interval
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	inline real integral_hom_2d(
		real(*f)(real, real),
		real a, real b,
		real c, real d,
		real f_max, PRNG& g,
		unsigned int N = 10000) {

		unsigned int N_inside = 0;

		for (unsigned int i = 0; i < N; ++i) {

			real x = rand_uniform(a, b, g);
			real y = rand_uniform(c, d, g);
			real z = rand_uniform(0, f_max, g);

			if(f(x, y) > z)
				N_inside++;
		}

		return (N_inside / static_cast<real>(N)) * (b - a) * (d - c) * f_max;
	}

}


#endif
