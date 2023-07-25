
/// @file rand_dist.h Random numbers following a pdf

#ifndef THEORETICA_RAND_DIST_H
#define THEORETICA_RAND_DIST_H


#include "../core/function.h"
#include "./prng.h"


namespace theoretica {


	/// A p.d.f sampling function taking as input
	/// the parameters of the distribution and
	/// a pseudorandom number generator.
	using pdf_sampling_function = real(*)(const std::vector<real>&, PRNG&);


	/// Generate a pseudorandom real number in [a, b] using a
	/// preexisting generator.
	/// @param a The lower extreme of the interval
	/// @param b The higher extreme of the interval
	/// @param g An already initialized pseudorandom number generator
	/// @param prec Precision parameters for the normalization, defaults
	/// to RAND_REAL_PREC.
	///
	/// The algorithm generates a random integer number, computes
	/// its modulus and divides it by prec:
	/// \f$x = \frac{(n mod p)}{2^p}\f$, where n is the random integer
	/// and p is the prec parameter
	inline real rand_uniform(real a, real b, PRNG& g, uint64_t prec = RAND_REAL_PREC) {

		// Generate a uniform random real number in [0, 1]
		real x = (g() % prec) / static_cast<real>(prec);

		// Transform to target interval
		return a + (b - a) * x;
	}


	/// Wrapper for rand_uniform(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_uniform(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("rand_uniform", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_uniform(theta[0], theta[1], g);
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
	/// MAX_TRYANDCATCH_ITER)
	/// @return A real number following the given pdf
	///
	/// Random real numbers are generated inside a rectangle
	/// defined by x1, x2, y1 and y2 following a uniform distribution.
	/// Only numbers below the pdf are returned.
	inline real rand_trycatch(stat_function f,
		const vec_buff& theta,
		real x1, real x2,
		real y1, real y2, PRNG& g,
		unsigned int max_iter = MAX_TRYANDCATCH_ITER) {

		real x;
		real y;

		unsigned int iter = 0;

		do {
			x = rand_uniform(x1, x2, g);
			y = rand_uniform(y1, y2, g);
			iter++;
		} while(y > f(x, theta) && iter <= max_iter);

		if(iter > max_iter) {
			TH_MATH_ERROR("rand_dist_tac", iter, NO_ALGO_CONVERGENCE);
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
	inline real rand_rejectsamp(
		stat_function f, const vec_buff& theta,
		real_function p, real_function Pinv,
		PRNG& g, unsigned int max_tries = 100) {

		for (unsigned int i = 0; i < max_tries; ++i) {

			// Generate a random number following
			// the p(x) probability distribution
			// by the inverse cumulative distribution function
			const real u_1 = rand_uniform(0, 1, g);
			const real x_p = Pinv(u_1);

			const real u_2 = rand_uniform(0, 1, g);

			// Accept the sample if f(x_p)/g(x_p) > u_2
			if(u_2 * p(x_p) < f(x_p, theta))
				return x_p;
		}

		TH_MATH_ERROR("rand_reject_sample", max_tries, NO_ALGO_CONVERGENCE);
		return nan();
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

			x = rand_uniform(-1, 1, g);
			y = rand_uniform(-1, 1, g);
			s = square(x) + square(y);

		} while(s >= 1 || s <= MACH_EPSILON);

		// Project the point
		s = sqrt(-2 * ln(s) / s);

		// Keep the second generated value for future calls
		spare = y * s;
		has_spare = true;

		return mean + sigma * x * s;
	}


	/// Generate a random number following a Gaussian distribution
	/// using the Box-Muller method.
	///
	/// @note This function may not be thread-safe as it uses
	/// static variables to keep track of spare generated values.
	inline real rand_gaussian_boxmuller(real mean, real sigma, PRNG& g) {

		static real spare;
		static bool has_spare = false;

		if(has_spare) {
			has_spare = false;
			return mean + spare * sigma;
		}

		// Generate a random point inside the unit circle
		
		const real x = rand_uniform(0, 1, g);
		const real y = rand_uniform(0, 1, g);

		const real x_transf = sqrt(-2 * ln(x));

		const real u = x_transf * cos(TAU * y);
		const real v = x_transf * sin(TAU * y);

		spare = v;
		has_spare = true;

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
	inline real rand_gaussian_clt(real mean, real sigma, PRNG& g) {

		// Fixed N = 12
		constexpr unsigned int N = 12;

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += rand_uniform(-1, 1, g);

		// f(u) = 1/2 (in [-1, 1])
		// E[u] = 0
		// sqrt(V[u]) = 1 / sqrt(3N) = 1 / 6

		return mean + (s / static_cast<real>(N)) * sigma * 6;
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
	inline real rand_gaussian_clt(
		real mean, real sigma,
		PRNG& g, unsigned int N) {

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += rand_uniform(-1, 1, g);

		// f(u) = 1/2 (in [-1, 1])
		// E[u] = 0
		// sqrt(V[u]) = 1 / sqrt(3N)

		return mean + (s / static_cast<real>(N)) * sigma * sqrt(3 * N);
	}


	/// Generate a random number following a Gaussian
	/// distribution using the best available algorithm.
	inline real rand_gaussian(real mean, real sigma, PRNG& g) {
		return rand_gaussian_polar(mean, sigma, g);
	}


	/// Wrapper for rand_gaussian(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_gaussian(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("rand_gaussian", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_gaussian(theta[0], theta[1], g);
	}


	/// Generate a random number following an exponential
	/// distribution using the quantile (inverse) function method.
	inline real rand_exponential(real lambda, PRNG& g) {

		if(abs(lambda) < MACH_EPSILON) {
			TH_MATH_ERROR("rand_exponential", lambda, DIV_BY_ZERO);
			return nan();
		}

		return -ln(1 - rand_uniform(0, 1, g)) / lambda;
	}


	/// Wrapper for rand_exponential(real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_exponential(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 1) {
			TH_MATH_ERROR("rand_exponential", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_exponential(theta[0], g);
	}


	/// Generate a random number following a Rayleigh
	/// distribution using the quantile (inverse) function method.
	inline real rand_rayleigh(real sigma, PRNG& g) {

		return sigma * sqrt(-2 * ln(1 - rand_uniform(0, 1, g)));
	}


	/// Wrapper for rand_rayleigh(real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_rayleigh(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 1) {
			TH_MATH_ERROR("rand_rayleigh", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_rayleigh(theta[0], g);
	}


	/// Generate a random number following a Cauchy
	/// distribution using the quantile (inverse) function method.
	inline real rand_cauchy(real mu, real alpha, PRNG& g) {

		return alpha * tan(PI * (rand_uniform(0, 1, g) - 0.5)) + mu;
	}


	/// Wrapper for rand_cauchy(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_cauchy(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("rand_cauchy", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_cauchy(theta[0], theta[1], g);
	}


	/// Generate a random number following a Laplace
	/// distribution using the quantile (inverse) function method.
	inline real rand_laplace(real mu, real b, PRNG& g) {

		real u = rand_uniform(0, 1, g);
		return mu - b * sgn(u) * ln(1 - 2 * abs(u));
	}


	/// Generate a random number following a Laplace
	/// distribution using the quantile (inverse) function method.
	inline real rand_laplace(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("rand_laplace", theta.size(), INVALID_ARGUMENT);
			return nan();
		}
		
		return rand_laplace(theta[0], theta[1], g);
	}


	/// Generate a random number following a Pareto
	/// distribution using the quantile (inverse) function method.
	inline real rand_pareto(real x_m, real alpha, PRNG& g) {

		return x_m / powf(1 - rand_uniform(0, 1, g), 1.0 / alpha);
	}


	/// Wrapper for rand_pareto(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	inline real rand_pareto(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("rand_pareto", theta.size(), INVALID_ARGUMENT);
			return nan();
		}

		return rand_pareto(theta[0], theta[1], g);
	}


	/// A probability density function sampler which
	/// generates pseudorandom numbers following
	/// asymptotically a given distribution.
	struct pdf_sampler {

		/// A p.d.f sampling function
		pdf_sampling_function f;

		/// The parameters of the target distribution
		std::vector<real> theta;

		/// A pseudorandom number generator
		PRNG& g;


		/// Initialize the sampler with the given parameters
		pdf_sampler(
			pdf_sampling_function f,
			const std::vector<real>& theta,
			PRNG& g) : f(f), theta(theta), g(g) {}


		/// Generate the next number
		inline real next() {
			return f(theta, g);
		}

		/// Generate the next number
		inline real operator()() {
			return next();
		}

		// Fill a vector with sampled points
		inline void fill(vec_buff& x, size_t N) {

			for (size_t i = 0; i < N; ++i)
				x[i] = next();
		}


		/// Returns a uniform distribution sampler
		static pdf_sampler uniform(real a, real b, PRNG& g) {
			return pdf_sampler(rand_uniform, {a, b}, g);
		}


		/// Returns a Gaussian distribution sampler
		static pdf_sampler gaussian(real mean, real sigma, PRNG& g) {
			return pdf_sampler(rand_gaussian, {mean, sigma}, g);
		}


		/// Returns an exponential distribution sampler
		static pdf_sampler exponential(real lambda, PRNG& g) {
			return pdf_sampler(rand_exponential, {lambda}, g);
		}


		/// Returns a Cauchy distribution sampler
		static pdf_sampler cauchy(real mu, real alpha, PRNG& g) {
			return pdf_sampler(rand_cauchy, {mu, alpha}, g);
		}


		/// Returns a Rayleigh distribution sampler
		static pdf_sampler rayleigh(real sigma, PRNG& g) {
			return pdf_sampler(rand_rayleigh, {sigma}, g);
		}


		/// Returns a Pareto distribution sampler
		static pdf_sampler pareto(real x_m, real alpha, PRNG& g) {
			return pdf_sampler(rand_pareto, {x_m, alpha}, g);
		}


		/// Returns a Laplace distribution sampler
		static pdf_sampler laplace(real mu, real b, PRNG& g) {
			return pdf_sampler(rand_laplace, {mu, b}, g);
		}

	};

	/// Metropolis algorithm for distribution sampling
	/// using a symmetric proposal distribution.
	///
	/// @param f The target distribution
	/// @param g A pdf_sampler already initialized to sample
	/// from the proposal distribution
	/// @param rnd An already initialized PRNG
	/// @param depth The number of iterations of the algorithm
	/// (defaults to METROPOLIS_DEPTH)
	inline real metropolis(
		real_function f, pdf_sampler& g,
		real x0, PRNG& rnd, unsigned int depth = METROPOLIS_DEPTH) {

		real current = x0, next;

		for(unsigned int i = 0; i < depth; i++) {
			
			// Computes the next step
			next = current + g();

			// Checks acceptance rate
			if(rand_uniform(0, 1, rnd) * f(current) <= f(next))
				current = next;
		}

		return current;
	}


	/// Metropolis algorithm for distribution sampling
	/// using a symmetric proposal distribution.
	/// This function uses the same PRNG as the proposal
	/// distribution sampler to generate uniform samples.
	///
	/// @param f The target distribution
	/// @param g A pdf_sampler already initialized to sample
	/// from the proposal distribution
	/// @param depth The number of iterations of the algorithm
	/// (defaults to METROPOLIS_DEPTH)
	inline real metropolis(real_function f, pdf_sampler& g,
		real x0, unsigned int depth = METROPOLIS_DEPTH) {
		return metropolis(f, g, x0, g.g, depth);
	}

}

#endif
