
/// @file sampling.h Sampling from probability distributions

#ifndef THEORETICA_SAMPLING_H
#define THEORETICA_SAMPLING_H

#include <tuple>

#include "../core/function.h"
#include "./random.h"


namespace theoretica {
namespace random {

	/// Generate a pseudorandom real number in [a, b] using a
	/// preexisting generator.
	///
	/// @param a The lower extreme of the interval
	/// @param b The higher extreme of the interval
	/// @param g An already initialized pseudorandom number generator
	/// @param prec Precision parameters for the normalization, defaults
	/// to PSEUDORANDOM_PREC.
	///
	/// The algorithm generates a random integer number, computes
	/// its modulus and divides it by prec:
	/// \f$x = \frac{(n mod p)}{2^p}\f$, where n is the random integer
	/// and p is the prec parameter
	template<typename PRNG>
	inline real uniform(PRNG& g, real a, real b, uint64_t prec = PSEUDORANDOM_PREC) {

		// Generate a uniform random real number in [0, 1]
		real x = (g() % prec) / static_cast<real>(prec);

		// Transform to target interval
		return a + (b - a) * x;
	}


	/// Fill an already initialized vector with uniformly distributed random numbers in [a, b].
	///
	/// @param x The vector to fill with random numbers
	/// @param a The lower extreme of the interval
	/// @param b The higher extreme of the interval
	/// @param g An already initialized pseudorandom number generator
	template<typename Vector, typename PRNG>
	inline void uniform(PRNG& g, Vector& x, real a, real b) {
		
		for (auto& x_i : x)
			x_i = uniform(g, a, b);
	}


	/// Generate a pseudorandom value following any
	/// probability distribution function using the
	/// Try-and-Catch (rejection) algorithm.
	///
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
	/// STATISTICS_TRYANDCATCH_ITER)
	/// @return A real number following the given pdf
	///
	/// Random real numbers are generated inside a rectangle
	/// defined by x1, x2, y1 and y2 following a uniform distribution.
	/// Only numbers below the pdf are returned.
	template<typename PRNG>
	inline real trycatch(
		stat_function f,
		const vec<real>& theta,
		real x1, real x2,
		real y1, real y2, PRNG& g,
		unsigned int max_iter = STATISTICS_TRYANDCATCH_ITER
	) {

		real x;
		real y;

		unsigned int iter = 0;

		do {
			x = random::uniform(g, x1, x2);
			y = random::uniform(g, y1, y2);
			iter++;
		} while(y > f(x, theta) && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("random::trycatch", iter, MathError::NoConvergence);
			return nan();
		}

		return x;
	}


	/// Generate a random number following any given distribution
	/// using rejection sampling.
	///
	/// @param f Target distribution
	/// @param theta The parameters of the target distribution
	/// @param p Proposal distribution
	/// @param Pinv Inverse cumulative function of the proposal distribution
	/// @param g An already initialized PRNG
	/// @param max_tries Maximum number of tries before
	/// stopping execution.
	template<typename PRNG>
	inline real rejectsamp(
		stat_function f, const vec<real>& theta,
		real_function p, real_function Pinv,
		PRNG& g, unsigned int max_tries = 100) {

		for (unsigned int i = 0; i < max_tries; ++i) {

			// Generate a random number following
			// the p(x) probability distribution
			// by the inverse cumulative distribution function
			const real u_1 = random::uniform(g, 0, 1);
			const real x_p = Pinv(u_1);

			const real u_2 = random::uniform(g, 0, 1); 

			// Accept the sample if f(x_p)/g(x_p) > u_2
			if(u_2 * p(x_p) < f(x_p, theta))
				return x_p;
		}

		TH_MATH_ERROR("random::rejectsamp", max_tries, MathError::NoConvergence);
		return nan();
	}


	/// Generate a random number following a Gaussian distribution
	/// using Marsaglia's polar method.
	template<typename PRNG>
	inline real gaussian_polar(PRNG& g, real mean, real sigma) {

		real x, y, s;

		// Generate a random point inside the unit circle
		do {

			x = random::uniform(g, -1, 1);
			y = random::uniform(g, -1, 1);
			s = square(x) + square(y);

		} while(s >= 1 || s <= MACH_EPSILON);

		// Project the point
		s = sqrt(-2 * ln(s) / s);

		return mean + sigma * x * s;
	}


	/// Generate a random number following a Gaussian distribution
	/// using the Box-Muller method.
	template<typename PRNG>
	inline real gaussian_boxmuller(PRNG& g, real mean, real sigma) {

		// Generate a random point inside the unit circle
		const real x = random::uniform(g, 0, 1);
		const real y = random::uniform(g, 0, 1);

		const real x_transf = sqrt(-2 * ln(x));
		const real u = x_transf * cos(TAU * y);

		return mean + sigma * u;
	}


	/// Generate a random number in a range
	/// following a Gaussian distribution by
	/// exploiting the Central Limit Theorem.
	/// @param mean The mean of the target distribution
	/// @param sigma The sigma of the target distribution
	/// @param g An already initialized PRNG
	///
	/// Exactly 12 real numbers in a range are generated
	/// and the mean is computed to get a single
	/// real number following (asymptotically) a
	/// Gaussian distribution.
	template<typename PRNG>
	inline real gaussian_clt(PRNG& g, real mean, real sigma) {

		// Fixed N = 12
		constexpr unsigned int N = 12;

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += random::uniform(g, -1, 1);

		// f(u) = 1/2 (in [-1, 1])
		// E[u] = 0
		// sqrt(V[u]) = 1 / sqrt(3N) = 1 / 6

		return mean + (s / N) * sigma * 6;
	}


	/// Generate a random number in a range
	/// following a Gaussian distribution by
	/// exploiting the Central Limit Theorem.
	/// @param mean The mean of the target distribution
	/// @param sigma The sigma of the target distribution
	/// @param g An already initialized PRNG
	/// @param N The number of random numbers to generate
	///
	/// Many real numbers in a range are generated
	/// and the mean is computed to get a single
	/// real number following (asymptotically) a
	/// Gaussian distribution.
	///
	/// @note This function uses a square root (th::sqrt)
	/// to rescale the output for variable N,
	/// the constant N implementation has better performance.
	template<typename PRNG>
	inline real gaussian_clt(
		PRNG& g, real mean, real sigma, unsigned int N
	) {

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += random::uniform(g, -1, 1);

		// f(u) = 1/2 (in [-1, 1])
		// E[u] = 0
		// sqrt(V[u]) = 1 / sqrt(3N)

		return mean + (s / N) * sigma * sqrt(3.0 * N);
	}


	/// Generate a random number following a Gaussian
	/// distribution using the best available algorithm.
	template<typename PRNG>
	inline real gaussian(PRNG& g, real mean, real sigma) {
		return random::gaussian_polar(g, mean, sigma);
	}


	/// Generate a random number following an exponential
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real exponential(PRNG& g, real lambda) {

		if(abs(lambda) < MACH_EPSILON) {
			TH_MATH_ERROR("random::exponential", lambda, MathError::DivByZero);
			return nan();
		}

		return -ln(1 - random::uniform(g, 0, 1)) / lambda;
	}


	/// Generate a random number following a Rayleigh
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real rayleigh(real sigma, PRNG& g) {

		return sigma * sqrt(-2 * ln(1 - random::uniform(g, 0, 1)));
	}


	/// Generate a random number following a Cauchy
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real cauchy(PRNG& g, real mu, real alpha) {

		return alpha * tan(PI * (random::uniform(g, 0, 1) - 0.5)) + mu;
	}


	/// Generate a random number following a Laplace
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real laplace(PRNG& g, real mu, real b) {

		const real u = random::uniform(g, 0, 0.5);
		return mu - b * rand_cointoss(g) * ln(1 - 2 * abs(u));
	}


	/// Generate a random number following a Pareto
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real pareto(PRNG& g, real x_m, real alpha) {

		return x_m / powf(1 - random::uniform(g, 0, 1), 1.0 / alpha);
	}


	/// Metropolis algorithm for distribution sampling
	/// using a symmetric proposal distribution.
	///
	/// @param pdf The target distribution
	/// @param proposal The proposal distribution, as an object
	/// satisfying the Distribution concept, such as stats::Gaussian.
	/// This object needs to provide the sample(PRNG&) method.
	/// @param rnd An already initialized PRNG
	/// @param depth The number of iterations of the algorithm
	/// (defaults to STATISTICS_METROPOLIS_DEPTH)
	template <
		typename PdfFunction,
		typename Distribution,
		typename ArgType,
		typename PRNG
	>
	inline ArgType metropolis(
		PdfFunction pdf, Distribution proposal,
		const ArgType& x0, PRNG& rnd,
		unsigned int depth = STATISTICS_METROPOLIS_DEPTH
	) {

		ArgType current = x0, next;

		for(unsigned int i = 0; i < depth; i++) {
			
			// Computes the next step
			next = current + proposal.sample(rnd);

			// Checks acceptance rate
			if(random::uniform(rnd, 0, 1) * pdf(current) <= pdf(next))
				current = next;
		}

		return current;
	}

}}

#endif
