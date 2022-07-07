
/// @file rand_dist.h Random numbers following a pdf

#ifndef THEORETICA_RAND_DIST_H
#define THEORETICA_RAND_DIST_H


#include "../core/function.h"
#include "./prng.h"


namespace theoretica {


	/// Generate a pseudorandom value following any
	/// probability distribution function using the
	/// Try-and-Catch (rejection) algorithm.
	/// @param f A probability distribution function
	/// @param theta The parameters of the pdf
	/// @param x1 The left extreme of the rectangle
	/// @param x2 The right extreme of the rectangle
	/// @param y1 The lower extreme of the rectangle
	/// @param y2 The upper extreme of the rectangle
	/// @param g An already initialized PRNG to use
	/// for number generation
	/// @param max_iter The maximum number of failed
	/// generations before stopping execution (defaults to
	/// MAX_TRYANDCATCH_ITER)
	/// @return A real number following the given pdf
	///
	/// Random real numbers are generated inside a rectangle
	/// defined by x1, x2, y1 and y2 following a uniform distribution.
	/// Only numbers below the pdf are returned.
	real rand_dist_tac(stat_function f,
		const vec_buff& theta,
		real x1, real x2,
		real y1, real y2, PRNG& g,
		unsigned int max_iter = MAX_TRYANDCATCH_ITER) {

		real x;
		real y;

		unsigned int iter = 0;

		do {
			x = rand_real(x1, x2, g);
			y = rand_real(y1, y2, g);
			iter++;
		} while(y > f(x, theta) && iter <= max_iter);

		if(iter > max_iter) {
			UMATH_ERROR("rand_dist_tac", iter, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Generate a random number following a Gaussian distribution
	/// using Marsaglia's polar method.
	///
	/// @note This function may not be thread-safe as it uses
	/// static variables to keep track of spare generated values.
	inline real rand_gaussian_polar(real mean, real sigma, PRNG& g) {

		static real spare;
		static bool has_spare = false;

		if(has_spare) {
			has_spare = false;
			return mean + spare * sigma;
		}

		real x, y, s;

		// Generate a random point inside the unit circle
		do {

			x = rand_real(-1, 1, g);
			y = rand_real(-1, 1, g);
			s = square(x) + square(y);

		} while(s >= 1 || s <= MACH_EPSILON);

		// Project the point onto the unit circumference
		s = sqrt(-2 * ln(s) / s);

		spare = y * s;
		has_spare = true;

		return mean + sigma * x * s;
	}


	/// Generate a random number in a range
	/// following a Gaussian distribution by
	/// exploiting the Central Limit Theorem.
	/// @param x1 Lower extreme of the interval of generation
	/// @param x2 Upper extreme of the interval of generation
	/// @param g An already initialized PRNG to use
	/// @param iter The number of random number to generate
	///
	/// Many real numbers in a range are generated
	/// and the mean is computed to get a single
	/// real number following (asymptotically) a
	/// Gaussian distribution.
	inline real rand_gaussian_clt(real x1, real x2, PRNG& g, unsigned int iter = 10) {

		real s = 0;
		for (unsigned int i = 0; i < iter; ++i)
			s += rand_real(x1, x2, g);

		return s / static_cast<real>(iter);
	}


	/// Generate a random number following a Gaussian
	/// distribution using the best available algorithm.
	inline real rand_gaussian(real mean, real sigma, PRNG& g) {
		return rand_gaussian_polar(mean, sigma, g);
	}


	/// Generate a random number following an exponential
	/// distribution using the quantile (inverse) function method.
	inline real rand_exponential(real lambda, PRNG& g) {

		if(abs(lambda) < MACH_EPSILON) {
			UMATH_ERROR("rand_exponential", lambda, DIV_BY_ZERO);
			return nan();
		}

		return -ln(1 - rand_real(0, 1, g)) / lambda;
	}


	/// Generate a random number following a Cauchy
	/// distribution using the quantile (inverse) function method.
	inline real rand_cauchy(real mu, real alpha, PRNG& g) {

		return alpha * tan(PI * (rand_real(0, 1, g) - 0.5)) + mu;
	}


	/// Generate a random number following a Pareto
	/// distribution using the quantile (inverse) function method.
	inline real rand_pareto(real x_m, real alpha, PRNG& g) {

		return x_m / powf(1 - rand_real(0, 1, g), 1.0 / alpha);
	}

}

#endif
