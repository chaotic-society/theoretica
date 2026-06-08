
/// @file sampling.h Sampling from probability distributions

#ifndef THEORETICA_SAMPLING_H
#define THEORETICA_SAMPLING_H

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
	inline real uniform(real a, real b, PRNG& g, uint64_t prec = PSEUDORANDOM_PREC) {

		// Generate a uniform random real number in [0, 1]
		real x = (g() % prec) / static_cast<real>(prec);

		// Transform to target interval
		return a + (b - a) * x;
	}


	/// Wrapper for random::uniform(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real uniform(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("random::uniform", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return uniform(theta[0], theta[1], g);
	}


	/// Fill an already initialized vector with uniformly distributed random numbers in [a, b].
	///
	/// @param x The vector to fill with random numbers
	/// @param a The lower extreme of the interval
	/// @param b The higher extreme of the interval
	/// @param g An already initialized pseudorandom number generator
	template<typename Vector, typename PRNG>
	inline void uniform(Vector& x, real a, real b, PRNG& g) {
		
		for (auto& x_i : x)
			x_i = uniform(a, b, g);
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
			x = random::uniform(x1, x2, g);
			y = random::uniform(y1, y2, g);
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
			const real u_1 = random::uniform(0, 1, g);
			const real x_p = Pinv(u_1);

			const real u_2 = random::uniform(0, 1, g);

			// Accept the sample if f(x_p)/g(x_p) > u_2
			if(u_2 * p(x_p) < f(x_p, theta))
				return x_p;
		}

		TH_MATH_ERROR("random::rejectsamp", max_tries, MathError::NoConvergence);
		return nan();
	}


	/// Generate a random number following a Gaussian distribution
	/// using Marsaglia's polar method.
	///
	/// @note This function may not be thread-safe as it uses
	/// static variables to keep track of spare generated values.
	template<typename PRNG>
	inline real gaussian_polar(real mean, real sigma, PRNG& g) {

		static real spare;
		static bool has_spare = false;

		if(has_spare) {
			has_spare = false;
			return mean + spare * sigma;
		}

		real x, y, s;

		// Generate a random point inside the unit circle
		do {

			x = random::uniform(-1, 1, g);
			y = random::uniform(-1, 1, g);
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
	template<typename PRNG>
	inline real gaussian_boxmuller(real mean, real sigma, PRNG& g) {

		static real spare;
		static bool has_spare = false;

		if(has_spare) {
			has_spare = false;
			return mean + spare * sigma;
		}

		// Generate a random point inside the unit circle
		
		const real x = random::uniform(0, 1, g);
		const real y = random::uniform(0, 1, g);

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
	template<typename PRNG>
	inline real gaussian_clt(real mean, real sigma, PRNG& g) {

		// Fixed N = 12
		constexpr unsigned int N = 12;

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += random::uniform(-1, 1, g);

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
	template<typename PRNG>
	inline real gaussian_clt(
		real mean, real sigma,
		PRNG& g, unsigned int N
	) {

		real s = 0;
		for (unsigned int i = 0; i < N; ++i)
			s += random::uniform(-1, 1, g);

		// f(u) = 1/2 (in [-1, 1])
		// E[u] = 0
		// sqrt(V[u]) = 1 / sqrt(3N)

		return mean + (s / static_cast<real>(N)) * sigma * sqrt(3 * N);
	}


	/// Generate a random number following a Gaussian
	/// distribution using the best available algorithm.
	template<typename PRNG>
	inline real gaussian(real mean, real sigma, PRNG& g) {
		return random::gaussian_polar(mean, sigma, g);
	}


	/// Wrapper for random::gaussian(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real gaussian(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("random::gaussian", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return random::gaussian(theta[0], theta[1], g);
	}


	/// Generate a random number following an exponential
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real exponential(real lambda, PRNG& g) {

		if(abs(lambda) < MACH_EPSILON) {
			TH_MATH_ERROR("random::exponential", lambda, MathError::DivByZero);
			return nan();
		}

		return -ln(1 - random::uniform(0, 1, g)) / lambda;
	}


	/// Wrapper for random::exponential(real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real exponential(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 1) {
			TH_MATH_ERROR("random::exponential", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return random::exponential(theta[0], g);
	}


	/// Generate a random number following a Rayleigh
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real rayleigh(real sigma, PRNG& g) {

		return sigma * sqrt(-2 * ln(1 - random::uniform(0, 1, g)));
	}


	/// Wrapper for random::rayleigh(real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real rayleigh(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 1) {
			TH_MATH_ERROR("random::rayleigh", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return random::rayleigh(theta[0], g);
	}


	/// Generate a random number following a Cauchy
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real cauchy(real mu, real alpha, PRNG& g) {

		return alpha * tan(PI * (random::uniform(0, 1, g) - 0.5)) + mu;
	}


	/// Wrapper for random::cauchy(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real cauchy(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("random::cauchy", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return random::cauchy(theta[0], theta[1], g);
	}


	/// Generate a random number following a Laplace
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real laplace(real mu, real b, PRNG& g) {

		const real u = random::uniform(0, 0.5, g);
		return mu - b * rand_cointoss(g) * ln(1 - 2 * abs(u));
	}


	/// Generate a random number following a Laplace
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real laplace(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("random::laplace", theta.size(), MathError::InvalidArgument);
			return nan();
		}
		
		return random::laplace(theta[0], theta[1], g);
	}


	/// Generate a random number following a Pareto
	/// distribution using the quantile (inverse) function method.
	template<typename PRNG>
	inline real pareto(real x_m, real alpha, PRNG& g) {

		return x_m / powf(1 - random::uniform(0, 1, g), 1.0 / alpha);
	}


	/// Wrapper for random::pareto(real, real, PRNG)
	///
	/// @param theta The parameters of the distribution
	/// @param g An already initialized PRNG
	template<typename PRNG>
	inline real pareto(const std::vector<real>& theta, PRNG& g) {

		if(theta.size() != 2) {
			TH_MATH_ERROR("random::pareto", theta.size(), MathError::InvalidArgument);
			return nan();
		}

		return random::pareto(theta[0], theta[1], g);
	}


	/// A probability density function sampler which
	/// generates pseudorandom numbers following
	/// asymptotically a given distribution
	/// \f$f(x; \vec \theta)\f$.
	template<typename PRNG>
	struct PdfSampler {

		/// A p.d.f sampling function
		real(*pdf)(const std::vector<real>&, PRNG&);

		/// The parameters of the target distribution
		std::vector<real> theta;

		/// A pseudorandom number generator
		PRNG& generator;


		/// Initialize the sampler with the given parameters
		PdfSampler(
			real(*pdf)(const std::vector<real>&, PRNG&),
			const std::vector<real>& theta,
			PRNG& generator) : pdf(pdf), theta(theta), generator(generator) {}


		/// Generate the next number
		inline real next() {
			return pdf(theta, generator);
		}

		/// Generate the next number
		inline real operator()() {
			return next();
		}

		/// Fill a vector with sampled points
		template<typename Vector>
		inline void fill(Vector& x, size_t N) {

			// Make sure that the buffer has enough space
			if (x.size() < N) {

				x.resize(N);

				if (x.size() < N) {
					TH_MATH_ERROR("PdfSampler::fill", N, MathError::InvalidArgument);
					algebra::vec_error(x);
					return;
				}
			}

			for (size_t i = 0; i < N; ++i)
				x[i] = next();
		}


		/// Fill a vector with sampled points
		template<typename Vector>
		inline void fill(Vector& x) {

			for (auto& x_i : x)
				x_i = next();
		}


		/// Stream the next generated number
		inline PdfSampler& operator>>(real& x) {
			x = next();
			return *this;
		}


		/// Returns a uniform distribution sampler
		static PdfSampler uniform(real a, real b, PRNG& generator) {
			return PdfSampler(random::uniform, {a, b}, generator);
		}


		/// Returns a Gaussian distribution sampler
		static PdfSampler gaussian(real mean, real sigma, PRNG& generator) {
			return PdfSampler(random::gaussian, {mean, sigma}, generator);
		}


		/// Returns an exponential distribution sampler
		static PdfSampler exponential(real lambda, PRNG& generator) {
			return PdfSampler(random::exponential, {lambda}, generator);
		}


		/// Returns a Cauchy distribution sampler
		static PdfSampler cauchy(real mu, real alpha, PRNG& generator) {
			return PdfSampler(random::cauchy, {mu, alpha}, generator);
		}


		/// Returns a Rayleigh distribution sampler
		static PdfSampler rayleigh(real sigma, PRNG& generator) {
			return PdfSampler(random::rayleigh, {sigma}, generator);
		}


		/// Returns a Pareto distribution sampler
		static PdfSampler pareto(real x_m, real alpha, PRNG& generator) {
			return PdfSampler(random::pareto, {x_m, alpha}, generator);
		}


		/// Returns a Laplace distribution sampler
		static PdfSampler laplace(real mu, real b, PRNG& generator) {
			return PdfSampler(random::laplace, {mu, b}, generator);
		}

	};

	/// Metropolis algorithm for distribution sampling
	/// using a symmetric proposal distribution.
	///
	/// @param pdf The target distribution
	/// @param g A PdfSampler already initialized to sample
	/// from the proposal distribution
	/// @param rnd An already initialized PRNG
	/// @param depth The number of iterations of the algorithm
	/// (defaults to STATISTICS_METROPOLIS_DEPTH)
	template <
		typename RealFunction,
		typename PRNG1,
		typename PRNG2
	>
	inline real metropolis(
		RealFunction pdf, PdfSampler<PRNG1>& g,
		real x0, PRNG2& rnd, unsigned int depth = STATISTICS_METROPOLIS_DEPTH) {

		real current = x0, next;

		for(unsigned int i = 0; i < depth; i++) {
			
			// Computes the next step
			next = current + g();

			// Checks acceptance rate
			if(random::uniform(0, 1, rnd) * pdf(current) <= pdf(next))
				current = next;
		}

		return current;
	}


	/// Metropolis algorithm for distribution sampling
	/// using a symmetric proposal distribution.
	/// This function uses the same PRNG as the proposal
	/// distribution sampler to generate uniform samples.
	///
	/// @param pdf The target distribution
	/// @param g A PdfSampler already initialized to sample
	/// from the proposal distribution
	/// @param depth The number of iterations of the algorithm
	/// (defaults to STATISTICS_METROPOLIS_DEPTH)
	template <
		typename RealFunction,
		typename PRNG
	>
	inline real metropolis(RealFunction pdf, PdfSampler<PRNG>& g,
		real x0, unsigned int depth = STATISTICS_METROPOLIS_DEPTH) {

		return metropolis(pdf, g, x0, g.generator, depth);
	}

}}

#endif
