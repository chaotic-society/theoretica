
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


namespace theoretica {


	namespace stats {

		/// Compute the likelihood of a distribution <f> with the given
		/// parameters <theta> and measures <X>
		/// @param X The dataset of the sample
		/// @param theta The parameters of the distribution
		/// @param f The statistical distribution function
		/// @result The likelihood of the given sample
		inline real likelihood(const vec<real>& X, const vec<real>& theta, stat_function f) {

			real res = 1;

			for (unsigned int i = 0; i < X.size(); ++i)
				res *= f(X[i], theta);

			return res;
		}


		/// Compute the log likelihood of a distribution <f> with the given
		/// parameters <theta> and measures <X>
		/// @param X The dataset of the sample
		/// @param theta The parameters of the distribution
		/// @param f The statistical distribution function
		/// @result The log-likelihood of the given sample
		inline real log_likelihood(const vec<real>& X, const vec<real>& theta, stat_function f) {

			real res = 0;

			for (unsigned int i = 0; i < X.size(); ++i)
				res += ln(f(X[i], theta));

			return res;
		}
	}


	/// @namespace theoretica::distribution Probability distribution functions
	namespace distribution {
	

		/// Gaussian distribution function
		inline real gaussian(real x, real X, real sigma) {

			return (1.0 / (sigma * SQRT2 * SQRTPI))
				* exp(-square(x - X) / (2 * square(sigma)));
		}


		/// Wrapper for distribution::gaussian(real, real, real)
		inline real gaussian(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::gaussian", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return gaussian(x, theta[0], theta[1]);
		}


		/// Bernoulli distribution
		inline real bernoulli(unsigned int k, real p) {

			return pow(p, k) * pow(1 - p, 1 - k);
		}


		/// Wrapper for distribution::bernoulli(unsigned int, real)
		inline real bernoulli(real k, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::bernoulli", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return bernoulli(static_cast<unsigned int>(k), theta[0]);

		}


		/// Poisson distribution
		inline real poisson(unsigned int k, real lambda) {

			return exp(-lambda) * pow(lambda, k) / (real) fact(k);
		}


		/// Wrapper for distribution::poisson(unsigned int, real)
		inline real poisson(real k, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::poisson", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return poisson(static_cast<unsigned int>(k), theta[0]);
		}


		/// Binomial distribution function
		inline real binomial(unsigned int nu, unsigned int n, real p) {

			return binomial_coeff(n, nu) *
				pow(p, nu) * pow(1 - p, n - nu);
		}


		/// Wrapper for distribution::binomial(unsigned int, unsigned int, real)
		inline real binomial(real nu, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::binomial", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return binomial(static_cast<unsigned int>(nu),
				static_cast<unsigned int>(theta[0]), theta[1]);
		}


		/// Multinomial distribution
		inline real multinomial(
			std::vector<unsigned int> x,
			unsigned int n,
			unsigned int k,
			std::vector<real> p) {

			if(x.size() != p.size() || x.size() != k) {
				TH_MATH_ERROR("distribution::multinomial", x.size(), INVALID_ARGUMENT);
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
		inline real chi_squared(real x, unsigned int k) {

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
		inline real chi_squared(real x, unsigned int k, real half_gamma_k) {

			if(x < 0)
				return 0;

			if(k % 2 == 0)
				return pow(x, int(k / 2) - 1) * exp(-x / 2.0)
					/ (pow(SQRT2, k) * half_gamma_k);
			else
				return pow(sqrt(x), k - 2) * exp(-x / 2.0)
					/ (pow(SQRT2, k) * half_gamma_k);
		}


		/// Wrapper for distribution::chi_squared(real, unsigned int)
		inline real chi_squared(real x, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::chi_squared", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return chi_squared(x, static_cast<unsigned int>(theta[0]));
		}


		/// Gamma distribution density function with parameters
		/// alpha (shape) and beta (rate)
		inline real gamma(real x, real alpha, real beta) {

			return powf(beta, alpha) * powf(x, alpha - 1)
				* exp(-beta * x) / special::gamma(alpha);
		}


		/// Wrapper for distribution::gamma(real, real, real)
		inline real gamma(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::gamma", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return gamma(x, theta[0], theta[1]);
		}


		/// Beta distribution density function with parameters
		/// alpha (shape) and beta (shape)
		inline real beta(real x, real alpha, real beta) {

			return powf(x, alpha - 1) * powf(1 - x, beta - 1)
				/ special::beta(alpha, beta);
		}


		/// Wrapper for distribution::beta(real, real, real)
		inline real beta(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::beta", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return beta(x, theta[0], theta[1]);
		}


		/// Student's t distribution
		inline real student(real x, unsigned int nu) {

			const real a = 1 + (x * x / nu);

			return (special::half_gamma(nu + 1) / special::half_gamma(nu))
					* pow(sqrt(a), -nu - 1) / (sqrt(nu) * SQRTPI);
		}


		/// Wrapper for distribution::student(real, unsigned int)
		inline real student(real x, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::student", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return student(x, static_cast<unsigned int>(theta[0]));
		}


		/// Log-normal distribution
		inline real log_normal(real x, real mu, real sigma) {

			return 1.0 / (2.50662827463 * sigma * x) *
				exp(-square(ln(x) - mu) / (2 * square(sigma)));
		}


		/// Wrapper for distribution::log_normal(real, real, real)
		inline real log_normal(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::log_normal", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return log_normal(x, theta[0], theta[1]);
		}


		/// Exponential distribution
		inline real exponential(real x, real lambda) {

			if(x < 0)
				return 0;

			return lambda * exp(-lambda * x);
		}


		/// Wrapper for distribution::exponential(real, real)
		inline real exponential(real x, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::exponential", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return exponential(x, theta[0]);
		}


		/// Rayleigh distribution
		inline real rayleigh(real x, real sigma) {

			if(x < 0)
				return 0;

			if(sigma < MACH_EPSILON) {
				TH_MATH_ERROR("distribution::rayleigh", sigma, DIV_BY_ZERO);
				return nan();
			}

			return x * exp(-square(x / sigma) / 2) / square(sigma);
		}


		/// Wrapper for distribution::rayleigh(real, real)
		inline real rayleigh(real x, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::rayleigh", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return rayleigh(x, theta[0]);
		}


		/// Cauchy distribution
		inline real cauchy(real x, real mu, real alpha) {

			return 1.0 / (PI * alpha * (1 + (square(x - mu) / square(alpha))));
		}


		/// Wrapper for distribution::cauchy(real, real, real)
		inline real cauchy(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::cauchy", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return cauchy(x, theta[0], theta[1]);
		}


		/// Breit Wigner distribution
		inline real breit_wigner(real x, real M, real Gamma) {

			return Gamma / (2 * PI * (square(x - M) + square(Gamma / 2.0)));
		}


		/// Wrapper for distribution::breit_wigner(real, real, real)
		inline real breit_wigner(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::breit_wigner", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return breit_wigner(x, theta[0], theta[1]);
		}


		/// Maxwell-Boltzmann distribution
		inline real maxwell(real x, real a) {

			return (SQRT2 / SQRTPI) * square(x) * exp(-square(x / a) / 2.0) / cube(a);
		}


		/// Wrapper for distribution::maxwell(real, real)
		inline real maxwell(real x, const vec<real>& theta) {

			if(theta.size() != 1) {
				TH_MATH_ERROR("distribution::maxwell", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return maxwell(x, theta[0]);
		}


		/// Laplace distribution
		inline real laplace(real x, real mu, real b) {

			return (1.0 / (2.0 * b)) * exp(-abs(x - mu) / b);
		}


		/// Laplace distribution
		inline real laplace(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::laplace", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return (1.0 / (2.0 * theta[1])) * exp(-abs(x - theta[0]) / theta[1]);
		}


		/// Pareto distribution
		inline real pareto(real x, real x_m, real alpha) {

			if(alpha <= 0) {
				TH_MATH_ERROR("distribution::pareto", alpha, OUT_OF_DOMAIN);
				return nan();
			}

			if(x < x_m)
				return 0;

			return alpha * powf(x_m, alpha) / powf(x, alpha + 1);
		}


		/// Wrapper for distribution::pareto(real, real, real)
		inline real pareto(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::pareto", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return pareto(x, theta[0], theta[1]);
		}


		/// Erlang distribution
		inline real erlang(real x, unsigned int k, real lambda) {

			if(!k) {
				TH_MATH_ERROR("distribution::erlang", k, INVALID_ARGUMENT);
				return nan();
			}

			return pow(lambda, k) * pow(x, k - 1)
					* exp(-lambda * x) / fact(k - 1);
		}


		/// Wrapper for distribution::erlang(real, unsigned int, real)
		inline real erlang(real x, const vec<real>& theta) {

			if(theta.size() != 2) {
				TH_MATH_ERROR("distribution::erlang", theta.size(), INVALID_ARGUMENT);
				return nan();
			}

			return erlang(x, static_cast<unsigned int>(theta[0]), theta[1]);
		}

	}
}

#endif
