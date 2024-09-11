
///
/// @file statistics.h Statistical functions
///

#ifndef THEORETICA_STATISTICS_H
#define THEORETICA_STATISTICS_H

#include "../core/constants.h"
#include "../algebra/vec.h"
#include "../core/real_analysis.h"
#include "../core/special.h"
#include "../calculus/integration.h"
#include "../calculus/gauss.h"


namespace theoretica {


	/// Compute the arithmetic mean of a set of values
	template<typename Dataset>
	inline real arithmetic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("arithmetic_mean", data.size(), DIV_BY_ZERO);
			return nan();
		}

		// Sum of x_i / N
		return sum(data) / (real) data.size();
	}


	/// Compute the arithmetic mean of a set of values
	/// Alias for arithmetic_mean
	template<typename Dataset>
	inline real mean(const Dataset& data) {
		return arithmetic_mean(data);
	}


	/// Compute the harmonic mean of a set of values
	template<typename Dataset>
	inline real harmonic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("harmonic_mean", data.size(), DIV_BY_ZERO);
			return nan();
		}

		real res = 0;

		for (unsigned int i = 0; i < data.size(); ++i) {

			if(data[i] == 0) {
				TH_MATH_ERROR("harmonic_mean", data[i], DIV_BY_ZERO);
				return nan();
			}

			res += 1.0 / data[i];
		}

		return static_cast<real>(data.size()) / res;
	}


	/// Compute the geometric mean of a set of values
	/// as \f$\sqrt[n]{\Pi_i x_i}\f$
	template<typename Dataset>
	inline real geometric_mean(const Dataset& data) {
		return root(product(data), data.size());
	}


	/// Compute the weighted mean of a set of values
	/// <data> and <weights> must have the same size
	template<typename Dataset1, typename Dataset2>
	inline real weighted_mean(const Dataset1& data, const Dataset2& weights) {

		// Sum of x_i * w_i / Sum of w_i
		return product_sum(data, weights) / sum(weights);
	}


	/// Compute the quadratic mean (Root Mean Square) of a set of values
	/// \f$m_q = \sqrt{x1^2 + x2^2 + ...}\f$
	template<typename Dataset>
	inline real quadratic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("quadratic_mean", data.size(), INVALID_ARGUMENT);
			return nan();
		}

		return sqrt(sum_squares(data) / data.size());
	}


	/// Compute the quadratic mean (Root Mean Square) of a set of values
	/// \f$m_q = \sqrt{x_1^2 + x_2^2 + ...}\f$
	/// @see quadratic_mean
	template<typename Dataset>
	inline real rms(const Dataset& data) {

		return quadratic_mean(data);
	}


	/// Computes the range of a data set
	/// defined as \f$x_{max} - {x_min}\f$
	template<typename Dataset>
	inline real range(const Dataset& data) {

		return max(data) - min(data);
	}


	/// Computes the maximum semidispersion of a data set
	/// defined as \f$(x_{max} - {x_min}) / 2\f$
	template<typename Dataset>
	inline real semidispersion(const Dataset& data) {

		return range(data) / 2.0;
	}


	/// Propagate error of a sum of values
	/// as sqrt(sigma_x^2 + sigma_y^2 + ...)
	template<typename Dataset>
	inline real propagate_sum(const Dataset& sigma) {

		return sqrt(sum_squares(sigma));
	}


	/// Propagate error of a product (or quotient) of values
	/// as \f$\sqrt{(sigma_x / x_{mean})^2 + (sigma_y / y_{mean})^2 + ...}\f$
	/// The result is the propagated relative error
	template<typename Dataset1, typename Dataset2>
	inline real propagate_product(const Dataset1& sigma, const Dataset2& mean) {

		if(sigma.size() != mean.size()) {
			TH_MATH_ERROR("propagate_product", sigma.size(), INVALID_ARGUMENT);
			return nan();
		}

		// Compute sum of squares of (i_sigma / i_mean)
		real s = 0;
		for (unsigned int i = 0; i < sigma.size(); ++i) {

			if(mean[i] == 0) {
				TH_MATH_ERROR("propagate_product", mean[i], DIV_BY_ZERO);
				return nan();
			}

			s += square(sigma[i] / abs(mean[i]));
		}

		return sqrt(s);
	}


	/// Total sum of squares (TSS) equal to \f$sum(square(x_i - x_{mean}))\f$,
	/// computed using Welford's one-pass method to improve numerical stability.
	template<typename Dataset>
	inline real total_sum_squares(const Dataset& X) {

		if(!X.size()) {
			TH_MATH_ERROR("total_sum_squares", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		// Running average
		real avg = X[0];

		// Total sum
		real s = 0.0;
		
		for (size_t i = 1; i < X.size(); ++i) {
			
			const real tmp = avg;

			avg = tmp + (X[i] - tmp) / (i + 1);
			s += (X[i] - tmp) * (X[i] - avg);
		}

		return s;
	}

	/// Total sum of squares (TSS)
	/// Computed as \f$\sum_i^n (x_i - x_{mean})^2 \f$
	/// @see total_sum_squares
	template<typename Dataset>
	inline real tss(const Dataset& X) {
		return total_sum_squares(X);
	}


	/// Compute the variance of a population
	template<typename Dataset>
	inline real variance(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("variance", data.size(), INVALID_ARGUMENT);
			return nan();
		}

		return total_sum_squares(data) / (real) data.size();
	}


	/// Compute the variance of a sample
	/// This function uses Bessel correction
	template<typename Dataset>
	inline real sample_variance(const Dataset& data) {

		if(data.size() < 2) {
			TH_MATH_ERROR("sample_variance", data.size(), INVALID_ARGUMENT);
			return nan();
		}

		// Bessel correction (N - 1)
		return total_sum_squares(data) / (real) (data.size() - 1);
	}


	/// Compute the standard deviation of a population
	template<typename Dataset>
	inline real standard_deviation(const Dataset& data) {
		return sqrt(variance(data));
	}


	/// Compute the standard deviation of a population
	template<typename Dataset>
	inline real stdev(const Dataset& data) {
		return standard_deviation(data);
	}


	/// Compute the standard deviation of a sample
	template<typename Dataset>
	inline real sample_standard_deviation(const Dataset& data) {
		return sqrt(sample_variance(data));
	}


	/// Compute the standard deviation of a sample
	template<typename Dataset>
	inline real smpl_stdev(const Dataset& data) {
		return sample_standard_deviation(data);
	}


	/// Compute the relative error on a population measure
	/// using standard deviation
	template<typename Dataset>
	inline real standard_relative_error(const Dataset& X) {

		real x_mean = mean(X);

		if(abs(x_mean) < MACH_EPSILON) {
			TH_MATH_ERROR("standard_relative_error", x_mean, DIV_BY_ZERO);
			return nan();
		}

		return standard_deviation(X) / abs(x_mean);
	}


	/// Compute the relative error on a sample measure
	/// using standard deviation
	template<typename Dataset>
	inline real sample_standard_relative_error(const Dataset& X) {

		real x_mean = mean(X);

		if(abs(x_mean) < MACH_EPSILON) {
			TH_MATH_ERROR("standard_relative_error", x_mean, DIV_BY_ZERO);
			return nan();
		}

		return sample_standard_deviation(X) / abs(x_mean);
	}


	/// Compute the standard deviation on the mean of a set of values
	template<typename Dataset>
	inline real mean_standard_deviation(const Dataset& data) {
		return sqrt(variance(data) / data.size());
	}


	/// Compute the standard deviation on the mean of a set of values
	template<typename Dataset>
	inline real stdom(const Dataset& data) {
		return mean_standard_deviation(data);
	}


	/// Compute the standard deviation on the mean of a set of measures
	/// Bessel correction is used in the calculation of the variance
	template<typename Dataset>
	inline real sample_mean_standard_deviation(const Dataset& data) {
		return sqrt(sample_variance(data) / data.size());
	}


	/// Compute the standard deviation on the mean of a set of measures
	/// Bessel correction is used in the calculation of the variance
	template<typename Dataset>
	inline real smpl_stdom(const Dataset& data) {
		return sample_mean_standard_deviation(data);
	}


	/// Compute the covariance of two sets of measures
	template<typename Dataset1, typename Dataset2>
	inline real covariance(const Dataset1& X, const Dataset2& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("covariance", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real s = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);

		for (unsigned int i = 0; i < X.size(); ++i)
			s += (X[i] - X_mean) * (Y[i] - Y_mean);

		return s / (real) X.size();
	}


	/// Compute the covariance between two sets of sample measures
	/// This function uses Bessel correction
	template<typename Dataset1, typename Dataset2>
	inline real sample_covariance(const Dataset1& X, const Dataset2& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("sample_covariance", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real s = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);

		for (unsigned int i = 0; i < X.size(); ++i)
			s += (X[i] - X_mean) * (Y[i] - Y_mean);

		// Bessel correction (N - 1)
		return s / (real) (X.size() - 1);
	}


	/// Pearson's correlation coefficient R for a population
	template<typename Dataset1, typename Dataset2>
	inline real correlation_coefficient(
		const Dataset1& X, const Dataset2& Y) {

		return covariance(X, Y) / (stdev(X) * stdev(Y));
	}


	/// Pearson's correlation coefficient r for a sample
	template<typename Dataset1, typename Dataset2>
	inline real sample_correlation_coefficient(
		const Dataset1& X, const Dataset2& Y) {

		return sample_covariance(X, Y) / (smpl_stdev(X) * smpl_stdev(Y));
	}


	/// Lag-n autocorrelation of a dataset
	/// @param X The dataset
	/// @param n The lag (defaults to lag-1)
	/// @return The lag-n autocorrelation of the given dataset
	template<typename Dataset>
	inline real autocorrelation(const Dataset& X, unsigned int n = 1) {

		if(X.size() < n) {
			TH_MATH_ERROR("autocorrelation", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		const real mu = mean(X);
		real num = 0;
		real den = square(X[0] - mu);

		for (unsigned int i = n; i < X.size(); ++i) {

			const real delta = X[i] - mu;
			num += delta * (X[i - n] - mu);
			den += square(delta);
		}

		return num / den;
	}


	/// Absolute deviation from the mean
	template<typename Dataset>
	inline real absolute_deviation(const Dataset& X) {

		real mu = mean(X);
		real res = 0;

		for (real x : X)
			res += abs(x - mu);

		return res / X.size();
	}


	/// Absolute deviation from the mean
	template<typename Dataset>
	inline real absdev(const Dataset& X) {
		return absolute_deviation(X);
	}


	/// Skewness of a dataset
	template<typename Dataset>
	inline real skewness(const Dataset& X) {

		const real mu = mean(X);
		const real sigma = smpl_stdev(X);
		real res = 0;

		for (real x : X)
			res += cube((x - mu) / sigma);

		return res / X.size();
	}


	/// Normalized Kurtosis of a dataset
	template<typename Dataset>
	inline real kurtosis(const Dataset& X) {

		const real mu = mean(X);
		const real sigma = smpl_stdev(X);
		real res = 0;

		for (real x : X)
			res += pow((x - mu) / sigma, 4);

		return (res / X.size()) - 3;
	}


	/// Compute the expectation of a given function with respect
	/// to a Gaussian distribution with the given parameters.
	/// This function uses Gauss-Hermite quadrature to compute
	/// the integral \f$\int_{-\infty}^{+\infty} g(x) e^{-x^2} dx\f$
	///
	/// @param mean The mean of the Gaussian distribution
	/// @param sigma The standard deviation of the Gaussian distribution
	/// @param g The function to compute the expectation of
	/// @return The Gaussian expectation of the given function
	template<typename RealFunction>
	inline real gaussian_expectation(RealFunction g, real mean, real sigma) {

		return integral_hermite(
			[=](real x) {
				return g(SQRT2 * sigma * x + mean);
			}
		) / SQRTPI;
	}


	/// Compute the Z-score of an observed value with
	/// respect to a Gaussian distribution with the
	/// given parameters
	///
	/// @param x The observed value
	/// @param mean The mean of the distribution
	/// @param sigma The standard deviation of the distribution
	/// @return The Z-score of x, computed as (x - mean) / sigma
	inline real z_score(real x, real mean, real sigma) {

		return (x - mean) / sigma;
	}


	/// Normalize a data set using Z-score normalization
	///
	/// @param The data set to normalize
	/// @return The normalized data set
	template<typename Dataset>
	inline Dataset normalize_z_score(const Dataset& X) {

		const real m = mean(X);
		const real s = sample_standard_deviation(X);

		return map([m, s](real x) { return z_score(x, m, s); }, X);
	}


	/// Compute the chi-square from the set of observed
	/// quantities, expected quantities and errors.
	/// The provided sets should all have the same size.
	///
	/// @param O The set of observed values
	/// @param E The set of expected values
	/// @param sigma The set of standard deviations on the observations
	/// @return The computed Chi-squared
	template<typename Dataset1, typename Dataset2, typename Dataset3>
	inline real chi_square(
		const Dataset1& O,
		const Dataset2& E,
		const Dataset3& sigma) {

		if(O.size() != E.size() || E.size() != sigma.size()) {
			TH_MATH_ERROR("chi_square", E.size(), INVALID_ARGUMENT);
			return nan();
		}

		real c_sqr = 0;

		for (unsigned int i = 0; i < O.size(); ++i) {

			if(abs(sigma[i]) < MACH_EPSILON) {
				TH_MATH_ERROR("chi_square", sigma[i], DIV_BY_ZERO);
				return nan();
			}

			c_sqr += square((O[i] - E[i]) / sigma[i]);
		}

		return c_sqr;
	}


	/// Compute the (right-tailed) p-value associated to a computed
	/// Chi-square value as the integral of the Chi-squared
	/// distribution from the given value to infinity (right-tailed).
	/// An equivalent integral is computed using Gauss-Laguerre quadrature: 
	/// \f$ p = \frac{e^{-X^2}}{2 \Gamma (k/2)} \int_0^{+\infty} (\sqrt{x + X^2})^{k - 2} e^{-x} dx \f$
	///
	/// @param chi_sqr The computed Chi-squared
	/// @param ndf Number of Degrees of Freedom
	/// @result The computed p-value
	/// @note The current implementation has reduced precision for 260 <= ndf < 1000
	/// because for ndf >= 260 the Gaussian approximation is used, which becomes
	/// more precise the higher the ndf.
	inline real pvalue_chi_squared(real chi_sqr, unsigned int ndf) {

		if(ndf == 0) {
			TH_MATH_ERROR("pvalue_chi_squared", ndf, INVALID_ARGUMENT);
			return nan();
		}

		// For ndf >= 260 use the Gaussian approximation
		// as the coefficients are not stable
		if(ndf >= 260) {
			
			const real new_x = (chi_sqr - ndf) / sqrt(2.0 * ndf);

			// For really low Chi-squared the Gaussian is
			// below tolerance value for integration
			if(new_x < 0) {

				if(new_x < -3)
					return 1 - integral_inf_riemann([=](real x) {
						return exp(-x * x / 2) / SQRTPI / SQRT2;
					}, -new_x, 1E-16, 25);

				return 0.5 + integral_romberg_tol([=](real x) {
					return exp(-x * x / 2) / SQRTPI / SQRT2;
				}, new_x, 0, 1E-16);
			} else {

				if(new_x > 3)
					return integral_inf_riemann([=](real x) {
						return exp(-x * x / 2) / SQRTPI / SQRT2;
					}, new_x, 1E-16, 25);

				return 0.5 - integral_romberg_tol([=](real x) {
					return exp(-x * x / 2) / SQRTPI / SQRT2;
				}, 0, new_x, 1E-16);
			}
		}

		// Compute the coefficient using a stable equivalent formula
		const real coeff = exp(-special::lngamma(ndf / 2.0) - chi_sqr / 2.0);

		// Use different methods when Gauss-Laguerre is not numerically stable
		if((ndf > 70 && chi_sqr < (ndf / 2.0))) {

			// Use equivalent formula around potential singularity
			real res = integral_romberg_tol([=](real x) {
				return pow(sqrt(x + chi_sqr / 2), ndf - 2) * exp(-x);
			}, 0, 1, 1E-12);

			res += integral_inf_riemann([=](real x) {
				return exp((ndf - 2) / 2.0 * ln(x + chi_sqr / 2) - x);
			}, 1, ndf / 2, 1E-12, 25);

			return coeff * res;
		}

		// Approximate the integral using Gauss-Laguerre quadrature
		return coeff * integral_gauss(
			[=](real x) {
				return pow(sqrt(x + chi_sqr / 2), ndf - 2);
		}, tables::laguerre_roots_16, tables::laguerre_weights_16, 16);
	}
  

	/// Compute the chi-square on a linearization
	template<typename Dataset1, typename Dataset2, typename Dataset3>
	inline real chi_square_linearization(
		const Dataset1& X, const Dataset2& Y,
		const Dataset3& sigma, real intercept, real slope) {

		if(X.size() != Y.size() || X.size() != sigma.size()) {
			TH_MATH_ERROR(
				"chi_square_linearization",
				X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real chi_squared = 0;
		for (unsigned int i = 0; i < X.size(); ++i) {

			if(abs(sigma[i]) <= MACH_EPSILON) {
				TH_MATH_ERROR("chi_square_linearization", sigma[i], DIV_BY_ZERO);
				return nan();
			}

			chi_squared += square((Y[i] - intercept - slope * X[i]) / sigma[i]);
		}

		return chi_squared;
	}


	/// Compute the reduced chi-squared on a linearization
	template<typename Dataset1, typename Dataset2, typename Dataset3>
	inline real reduced_chi_square_linearization(
		const Dataset1& X, const Dataset2& Y,
		const Dataset3& sigma, real intercept, real slope) {

		if(Y.size() <= 2) {
			TH_MATH_ERROR("reduced_chi_square_linearization",
				Y.size(), INVALID_ARGUMENT);
			return nan();
		}

		// Divide by degrees of freedom (N - 2)
		return chi_square_linearization(X, Y, sigma, intercept, slope)
			/ (real) (Y.size() - 2);
	}

}

#endif
