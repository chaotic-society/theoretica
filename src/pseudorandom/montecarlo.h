
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
		PRNG& g, unsigned int N = 1000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i)			
			sum_y += f(rand_uniform(a, b, g));

		return (b - a) * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Crude Monte Carlo integration
	/// @param f The function to integrate
	/// @param extremes A vector of the extremes of integration
	/// @param g An already initialized PRNG
	/// @param N The number of points to generate
	template<unsigned int S>
	inline real integral_crude(
		real(*f)(vec<S>), vec<S, vec2> extremes,
		PRNG& g, unsigned int N = 1000) {

		real sum_y = 0;

		// Sample the function at random points in the integration region
		for (unsigned int i = 0; i < N; ++i) {
		
			vec<S> v;
			for (unsigned int k = 0; k < S; ++k)
				v[k] = rand_uniform(extremes[k][0], extremes[k][1], g);

			sum_y += f(v);
		}

		// Volume of the integration region
		real volume = 1;
		for (unsigned int i = 0; i < S; ++i)
			volume *= abs(extremes[i][1] - extremes[i][0]);

		return volume * sum_y / static_cast<real>(N);
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
		unsigned int N = 1000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i)
			sum_y += f(a + qrand_weyl(i) * (b - a));

		return (b - a) * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Crude
	/// Quasi-Monte Carlo integration by sampling
	/// from the Weyl sequence
	/// @param f The function to integrate
	/// @param extremes A vector of the extremes of integration
	/// @param alpha A vector of the irrational numbers to use for the Weyl sequence
	/// @param N The number of points to generate
	template<unsigned int S>
	inline real integral_quasi_crude(
		real(*f)(vec<S>), vec<S, vec2> extremes,
		unsigned int N, vec<S> alpha) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i) {
		
			vec<S> v;
			for (unsigned int k = 0; k < S; ++k)
				v[k] = extremes[k][0] + (qrand_weyl(i, alpha[k]) * (extremes[k][1] - extremes[k][0]));

			sum_y += f(v);
		}

		real volume = 1;
		for (unsigned int i = 0; i < S; ++i)
			volume *= abs(extremes[i][1] - extremes[i][0]);

		return volume * sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Crude
	/// Quasi-Monte Carlo integration by sampling
	/// from the Weyl sequence
	/// @param f The function to integrate
	/// @param extremes A vector of the extremes of integration
	/// @param alpha An irrational number
	/// @param N The number of points to generate
	template<unsigned int S>
	inline real integral_quasi_crude(
		real(*f)(vec<S>), vec<S, vec2> extremes,
		unsigned int N = 1000, real alpha = 0) {

		if(alpha == 0) {

			polynomial<> P = polynomial<real>::monomial(1, S + 1);
			P[0] = -1;
			P[1] = -1;

			// Find the only positive root of the polynomial
			// x^s+1 - x - 1 = 0

			alpha = 1.0 / approx_root_bisection(
				[P](real x) {
					return P(x);
				}, 0, 2);
		}

		vec<S> v;
		for (unsigned int i = 0; i < S; ++i)
			v[i] = pow(alpha, i + 1);

		return integral_quasi_crude(f, extremes, N, v);
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
		PRNG& g, unsigned int N = 1000) {

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
		unsigned int N = 1000) {

		unsigned int N_inside = 0;

		for (unsigned int i = 0; i < N; ++i) {

			vec2 v = qrand_weyl2(i);
		
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
		unsigned int N = 1000) {

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


	/// Approximate an integral by using Crude Monte Carlo integration with
	/// importance sampling.
	/// @param f The function to integrate
	/// @param g The importance function (normalized)
	/// @param Ginv The inverse of the primitive of g, with domain [0, 1]
	/// @param gen An already initialized PRNG
	/// @param N The number of points to generate
	inline real integral_impsamp(
		real_function f, real_function g, real_function Ginv,
		PRNG& gen, unsigned int N = 1000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i) {
			const real z = Ginv(rand_uniform(0, 1, gen));
			sum_y += f(z) / g(z);
		}		

		return sum_y / static_cast<real>(N);
	}


	/// Approximate an integral by using Crude Quasi-Monte Carlo integration with
	/// importance sampling, using the Weyl sequence.
	/// @param f The function to integrate
	/// @param g The importance function (normalized)
	/// @param Ginv The inverse of the primitive of g, with domain [0, 1]
	/// @param gen An already initialized PRNG
	/// @param N The number of points to generate
	inline real integral_quasi_impsamp(
		real_function f, real_function g, real_function Ginv,
		unsigned int N = 1000) {

		real sum_y = 0;

		for (unsigned int i = 0; i < N; ++i) {
			const real z = Ginv(qrand_weyl(i + 1));
			sum_y += f(z) / g(z);
		}		

		return sum_y / static_cast<real>(N);
	}

}


#endif
