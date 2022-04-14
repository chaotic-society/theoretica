
/// @file rand_dist.h Random numbers following a pdf

#ifndef UROBORO_RAND_DIST_H
#define UROBORO_RAND_DIST_H


#include "../function.h"
#include "./prng.h"


namespace uroboro {


	/// Generate a pseudorandom value following any
	/// probability distribution function using the
	/// Try-and-Catch algorithm.
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
		for (int i = 0; i < iter; ++i)
			s += rand_real(x1, x2, g);

		return s / static_cast<real>(iter);
	}

}

#endif
