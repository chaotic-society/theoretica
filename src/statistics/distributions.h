
///
/// @file distributions.h Probability distribution functions
///

#ifndef THEORETICA_DISTRIBUTIONS
#define THEORETICA_DISTRIBUTIONS

#include "../core/real_analysis.h"
#include "../algebra/vec.h"
#include "./statistics.h"
#include "../core/function.h"
#include "../core/special.h"
#include "../random/sampling.h"


namespace theoretica {
namespace stats {

	/// Compute the likelihood of a probability distribution
	/// with the given parameters for a dataset of measures.
	/// 
	/// @param X The dataset of the sample
	/// @param f The statistical distribution function
	/// @param params The parameters of the distribution as a variadic template list.
	/// @result The likelihood of the given sample
	template <
		typename Vector,
		typename Function,
		typename ...Params
	>
	inline real likelihood(const Vector& X, Function f, Params... params) {

		real res = 1;

		for (unsigned int i = 0; i < X.size(); ++i)
			res *= f(X[i], params...);

		return res;
	}


	/// Compute the log-likelihood of a probability distribution
	/// with the given parameters for a dataset of measures.
	/// 
	/// @param X The dataset of the sample
	/// @param f The statistical distribution function
	/// @param params The parameters of the distribution as a variadic template list.
	/// @result The log likelihood of the given sample
	template <
		typename Vector,
		typename Function,
		typename ...Params
	>
	inline real log_likelihood(const Vector& X, Function f, Params... params) {

		real res = 0;

		for (unsigned int i = 0; i < X.size(); ++i)
			res += ln(f(X[i], params...));

		return res;
	}


	// Statistical distributions


	/// Uniform distribution function
	inline real uniform(real x, real a, real b) {

		if(x < a || x > b)
			return 0;

		return 1.0 / (b - a);
	}


	/// Gaussian distribution function
	inline real gaussian(real x, real X, real sigma) {

		return (1.0 / (sigma * SQRT2 * SQRTPI))
			* exp(-square(x - X) / (2 * square(sigma)));
	}


	/// Bernoulli distribution
	inline real bernoulli(unsigned int k, real p) {

		return pow(p, k) * pow(1 - p, 1 - k);
	}


	/// Poisson distribution
	inline real poisson(unsigned int k, real lambda) {

		return exp(-lambda) * pow(lambda, k) / (real) fact(k);
	}


	/// Binomial distribution function
	inline real binomial(unsigned int nu, unsigned int n, real p) {

		return binomial_coeff(n, nu) *
			pow(p, nu) * pow(1 - p, n - nu);
	}


	/// Multinomial distribution
	inline real multinomial(
		std::vector<unsigned int> x,
		unsigned int n,
		unsigned int k,
		std::vector<real> p) {

		if(x.size() != p.size() || x.size() != k) {
			TH_MATH_ERROR("stats::multinomial", x.size(), MathError::InvalidArgument);
			return nan();
		}

		real res = fact(n);

		for (unsigned int i = 0; i < k; ++i)
			res *= pow(p[i], x[i]) / fact(x[i]);

		return res;
	}


	/// Chi-squared distribution
	///
	/// @param x The point to evaluate the distribution at
	/// @param k The number of degrees of freedom
	inline real chi2(real x, unsigned int k) {

		if(x < 0)
			return 0;

		const real coeff = special::half_gamma(k);

		if(k % 2 == 0)
			return pow(x, int(k / 2) - 1) * exp(-x / 2.0)
				/ (pow(SQRT2, k) * coeff);
		else
			return pow(sqrt(x), k - 2) * exp(-x / 2.0)
				/ (pow(SQRT2, k) * coeff);
	}


	/// Chi-squared distribution
	///
	/// @param x The point to evaluate the distribution at
	/// @param k The number of degrees of freedom
	/// @param half_gamma_k A precomputed value of special::half_gamma(k)
	///
	/// This function accepts a precomputed value of special::half_gamma(k)
	/// for repeated evaluation of the distribution. You can compute it once
	/// and reuse the same result for constant ndf.
	inline real chi2(real x, unsigned int k, real half_gamma_k) {

		if(x < 0)
			return 0;

		if(k % 2 == 0)
			return pow(x, int(k / 2) - 1) * exp(-x / 2.0)
				/ (pow(SQRT2, k) * half_gamma_k);
		else
			return pow(sqrt(x), k - 2) * exp(-x / 2.0)
				/ (pow(SQRT2, k) * half_gamma_k);
	}


	/// Gamma distribution density function with parameters
	/// alpha (shape) and beta (rate)
	inline real gamma(real x, real alpha, real beta) {

		return powf(beta, alpha) * powf(x, alpha - 1)
			* exp(-beta * x) / special::gamma(alpha);
	}


	/// Beta distribution density function with parameters
	/// alpha (shape) and beta (shape)
	inline real beta(real x, real alpha, real beta) {

		return powf(x, alpha - 1) * powf(1 - x, beta - 1)
			/ special::beta(alpha, beta);
	}


	/// Student's t distribution
	inline real student(real x, unsigned int nu) {

		const real a = 1 + (x * x / nu);

		return (special::half_gamma(nu + 1) / special::half_gamma(nu))
			* pow(sqrt(a), -nu - 1) / (sqrt(nu) * SQRTPI);
	}


	/// Log-normal distribution
	inline real log_normal(real x, real mu, real sigma) {

		return 1.0 / (2.50662827463 * sigma * x) *
			exp(-square(ln(x) - mu) / (2 * square(sigma)));
	}


	/// Exponential distribution
	inline real exponential(real x, real lambda) {

		if(x < 0)
			return 0;

		return lambda * exp(-lambda * x);
	}


	/// Rayleigh distribution
	inline real rayleigh(real x, real sigma) {

		if(x < 0)
			return 0;

		if(sigma < MACH_EPSILON) {
			TH_MATH_ERROR("stats::rayleigh", sigma, MathError::DivByZero);
			return nan();
		}

		return x * exp(-square(x / sigma) / 2) / square(sigma);
	}


	/// Cauchy distribution
	inline real cauchy(real x, real mu, real alpha) {

		return 1.0 / (PI * alpha * (1 + (square(x - mu) / square(alpha))));
	}


	/// Breit Wigner distribution
	inline real breit_wigner(real x, real M, real Gamma) {

		return Gamma / (2 * PI * (square(x - M) + square(Gamma / 2.0)));
	}


	/// Maxwell-Boltzmann distribution
	inline real maxwell(real x, real a) {

		return (SQRT2 / SQRTPI) * square(x) * exp(-square(x / a) / 2.0) / cube(a);
	}


	/// Laplace distribution
	inline real laplace(real x, real mu, real b) {

		return (1.0 / (2.0 * b)) * exp(-abs(x - mu) / b);
	}


	/// Pareto distribution
	inline real pareto(real x, real x_m, real alpha) {

		if(alpha <= 0) {
			TH_MATH_ERROR("stats::pareto", alpha, MathError::OutOfDomain);
			return nan();
		}

		if(x < x_m)
			return 0;

		return alpha * powf(x_m, alpha) / powf(x, alpha + 1);
	}


	/// Erlang distribution
	inline real erlang(real x, unsigned int k, real lambda) {

		if(!k) {
			TH_MATH_ERROR("stats::erlang", k, MathError::InvalidArgument);
			return nan();
		}

		return pow(lambda, k) * pow(x, k - 1) * exp(-lambda * x) / fact(k - 1);
	}



	/// @class Gaussian
	/// A class representing the Gaussian distribution,
	/// with methods to evaluate the distribution and
	/// sample from it. This class can be used to represent
	/// a Gaussian distribution and fulfills the Distribution concept,
	/// providing a pdf() and a sample() method.
	class Gaussian {
	private:
		real mu {0.0};
		real sigma {1.0};

	public:

		/// Construct the standard Gaussian distribution
		/// with zero mean and standard deviation 1.
		Gaussian() = default;

		/// Construct a Gaussian distribution with the given moments.
		/// @param mu The mean of the distribution
		/// @param sigma The standard deviation of the distribution
		Gaussian(real mu, real sigma) : mu(mu), sigma(sigma) {}

		/// Evaluate the Gaussian distribution at a given point.
		real pdf(real x) const {
			return stats::gaussian(x, mu, sigma);
		}

		/// Evaluate the Gaussian distribution at a given point.
		/// Equivalent to eval.
		real operator()(real x) const {
			return pdf(x);
		}

		/// Sample a single value from the Gaussian distribution.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @return The sampled value, following the Gaussian distribution.
		template<typename PRNG>
		real sample(PRNG& g) const {
			return random::gaussian(g, mu, sigma);
		}

		/// Sample a vector from the Gaussian distribution,
		/// filling an already initialized vector.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @param v The vector to fill with sampled values.
		template<typename PRNG, typename Vector>
		Vector& sample(PRNG& g, Vector& v) const {
			
			for (auto& x : v)
				x = random::gaussian(g, mu, sigma);
			
			return v;
		}
	};


	/// @class Uniform
	/// A class representing the uniform distribution,
	/// with methods to evaluate the distribution and
	/// sample from it. This class can be used to represent
	/// a uniform distribution and fulfills the Distribution concept,
	/// providing a pdf() and a sample() method.
	class Uniform {
	private:
		real a {0.0};
		real b {1.0};


	public:

		/// Construct the standard uniform distribution
		/// with lower bound 0 and upper bound 1.
		Uniform() = default;

		/// Construct a uniform distribution with the given bounds.
		/// @param a The lower bound of the distribution
		/// @param b The upper bound of the distribution
		Uniform(real a, real b) : a(a), b(b) {}

		/// Evaluate the uniform distribution at a given point.
		real pdf(real x) const {
			return stats::uniform(x, a, b);
		}

		/// Evaluate the uniform distribution at a given point.
		/// Equivalent to eval.
		real operator()(real x) const {
			return pdf(x);
		}

		/// Sample a single value from the uniform distribution.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @return The sampled value, following the uniform distribution.
		template<typename PRNG>
		real sample(PRNG& g) const {
			return random::uniform(g, a, b);
		}

		/// Sample a vector from the uniform distribution,
		/// filling an already initialized vector.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @param v The vector to fill with sampled values.
		template<typename PRNG, typename Vector>
		Vector& sample(PRNG& g, Vector& v) const {
			
			for (auto& x : v)
				x = random::uniform(g, a, b);
			
			return v;
		}
	};


	/// @class Exponential
	/// A class representing the exponential distribution,
	/// with methods to evaluate the distribution and
	/// sample from it. This class can be used to represent
	/// a uniform distribution and fulfills the Distribution concept,
	/// providing a pdf() and a sample() method.
	class Exponential {
	private:
		real lambda {1.0};

	public:

		/// Construct the standard exponential distribution
		/// with rate parameter 1.
		Exponential() = default;

		/// Construct an exponential distribution with the given bounds.
		/// @param a The lower bound of the distribution
		/// @param b The upper bound of the distribution
		Exponential(real lambda) : lambda(lambda) {}

		/// Evaluate the exponential distribution at a given point.
		real pdf(real x) const {
			return stats::exponential(x, lambda);
		}

		/// Evaluate the exponential distribution at a given point.
		/// Equivalent to eval.
		real operator()(real x) const {
			return pdf(x);
		}

		/// Sample a single value from the exponential distribution.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @return The sampled value, following the exponential distribution.
		template<typename PRNG>
		real sample(PRNG& g) const {
			return random::exponential(g, lambda);
		}

		/// Sample a vector from the exponential distribution,
		/// filling an already initialized vector.
		/// @param g An already initialized pseudorandom number generator (e.g. XoshiroPrng).
		/// @param v The vector to fill with sampled values.
		template<typename PRNG, typename Vector>
		Vector& sample(PRNG& g, Vector& v) const {
			for (auto& x : v)
				x = random::exponential(g, lambda);
			return v;
		}
	};

}}

#endif
