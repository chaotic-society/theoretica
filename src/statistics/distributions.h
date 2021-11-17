#ifndef UROBORO_DISTRIBUTIONS
#define UROBORO_DISTRIBUTIONS

#include "../common.h"
#include "./statistics.h"


namespace uroboro {

	namespace distribution {
	

		// Gaussian distribution function
		inline real gaussian(real x, real X, real sigma) {

			return (1.0 / (sigma *
				uroboro::sqrt(2 * PI))) * uroboro::exp(-square(x - X) / (2 * square(sigma)));
		}


		// Gaussian distribution function calculated on a sample of measures
		inline real gaussian(real x, const vec_buff& data) {

			return gaussian(x, mean(data),
				sample_mean_standard_deviation(data));
		}


		// Bernoulli distribution
		inline real bernoulli(unsigned int k, real p) {

			return pow(p, k) * pow(1 - p, 1 - k);
		}


		// Poisson distribution
		inline real poisson(unsigned int k, real lambda) {

			return exp(-lambda) * pow(lambda, k) / (real) fact(k);
		}


		// Binomial distribution function
		inline real binomial(unsigned int nu, unsigned int n, real p) {

			return binomial_coeff(n, nu) *
				uroboro::pow(p, nu) * uroboro::pow(1 - p, n - nu);
		}


		// Log-normal distribution
		inline real log_normal(real x, real mu, real sigma) {

			return 1.0 / (2.50662827463 * sigma * x) *
				exp(-square(ln(x) - mu) / (2 * square(sigma)));
		}


		// Exponential distribution
		inline real exponential(real x, real lambda) {

			if(x < 0)
				return 0;

			return lambda * exp(-lambda * x);
		}


		// Cauchy distribution
		inline real cauchy(real x, real mu, real alfa) {

			return 1.0 / (PI * alfa * (1 + (square(x - mu) / square(alfa))));
		}


		// Breit Wigner distribution
		inline real breit_wigner(real x, real M, real Gamma) {

			return Gamma / (2 * PI * (square(x - M) + square(Gamma / 2.0)));
		}

	}
}

#endif
