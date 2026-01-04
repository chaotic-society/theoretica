
///
/// @file statistics.h Statistical functions
///

#ifndef THEORETICA_STATISTICS_H
#define THEORETICA_STATISTICS_H

#include "../core/constants.h"
#include "../core/real_analysis.h"
#include "../core/special.h"
#include "../calculus/integral.h"
#include "../calculus/gauss.h"
#include "../core/dataset.h"


namespace theoretica {


	/// @namespace theoretica::stats Statistical functions
	namespace stats {


		/// Compute the mean of a dataset
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The mean of the dataset
		template<typename Dataset>
		inline real mean(const Dataset& X) {
			return arithmetic_mean(X);
		}


		/// Computes the range of a data set, defined as \f$x_{max} - {x_min}\f$
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The range of the values of the dataset
		template<typename Dataset>
		inline real range(const Dataset& X) {

			return max(X) - min(X);
		}


		/// Computes the maximum semidispersion of a data set
		/// defined as \f$(x_{max} - {x_min}) / 2\f$.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The maximum semidispersion of the dataset
		template<typename Dataset>
		inline real semidispersion(const Dataset& X) {

			return range(X) / 2.0;
		}


		/// Propagate the error over a sum of random variables under quadrature, as
		/// \f$\sqrt{\sum_{i = 1}^n \sigma_i^2}\f$, where each \f$\sigma_i\f$
		/// corresponds to the standard deviation of a variable.
		/// The random variables are assumed to be statistically independent.
		///
		/// @param sigma The vector of standard deviations
		/// @return The propagated error over the sum
		template<typename Dataset>
		inline real propagate_sum(const Dataset& sigma) {

			return sqrt(sum_squares(sigma));
		}


		/// Propagate the error over a product of random variables under quadrature, as
		/// \f$\sqrt{\sum_{i = 1}} (\sigma_i / \mu_i)^2}\f$,
		/// where each \f$\sigma_i\f$ corresponds to the standard deviation of a variable.
		/// The random variables are assumed to be statistically independent and the
		/// result is the relative error over the product.
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @param sigma The vector of standard deviations
		/// @param mean The vector of the mean values
		/// @return The propagated relative error over the product
		template<typename Dataset1, typename Dataset2>
		inline real propagate_product(const Dataset1& sigma, const Dataset2& mean) {

			if(sigma.size() != mean.size()) {
				TH_MATH_ERROR("propagate_product", sigma.size(), MathError::InvalidArgument);
				return nan();
			}

			// Compute sum of squares of (i_sigma / i_mean)
			real s = 0;
			for (unsigned int i = 0; i < sigma.size(); ++i) {

				if(mean[i] == 0) {
					TH_MATH_ERROR("propagate_product", mean[i], MathError::DivByZero);
					return nan();
				}

				s += square(sigma[i] / abs(mean[i]));
			}

			return sqrt(s);
		}


		/// Compute the total sum of squares (TSS) of a given dataset
		/// as \f$sum(square(x_i - x_{mean}))\f$ using Welford's one-pass method.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset to compute the TSS of
		/// @return The total sum of squares of the given dataset
		template<typename Dataset>
		inline real total_sum_squares(const Dataset& X) {

			if(!X.size()) {
				TH_MATH_ERROR("total_sum_squares", X.size(), MathError::InvalidArgument);
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


		/// Compute the variance given a dataset and the number
		/// of constraints. Welford's one-pass method is used. The number
		/// of constraints defaults to 1, applying Bessel's correction.
		/// A value of 0 may be used to compute the population variance.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @param constraints The number of constraints, defaults to 1
		/// @return The variance of the dataset
		template<typename Dataset>
		inline real variance(const Dataset& X, unsigned int constraints = 1) {

			if(X.size() <= constraints) {
				TH_MATH_ERROR("variance", X.size(), MathError::InvalidArgument);
				return nan();
			}

			return total_sum_squares(X) / (X.size() - constraints);
		}


		/// Compute the mean and the variance of a dataset
		/// in a single pass, using Welford's method, with the given number of constraints
		/// (defaults to 1 for Bessel's correction).
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @param out_mean A reference to overwrite with the computed mean
		/// @param out_variance A reference to overwrite with the computed variance
		/// @param constraints The number of constraints (defaults to 1)
		template<typename Dataset>
		inline void moments2(
			const Dataset& X, real& out_mean,
			real& out_variance, unsigned int constraints = 1) {

			if(X.size() <= constraints) {
				TH_MATH_ERROR("total_sum_squares", X.size(), MathError::InvalidArgument);
				out_mean = nan();
				out_variance = nan();
				return;
			}

			// Running average
			real avg = X[0];

			// Total sum
			real tss = 0.0;
			
			for (size_t i = 1; i < X.size(); ++i) {
				
				const real tmp = avg;

				avg = tmp + (X[i] - tmp) / (i + 1);
				tss += (X[i] - tmp) * (X[i] - avg);
			}

			out_mean = avg;
			out_variance = tss / (X.size() - constraints);
		}


		/// Compute the standard deviation given a dataset and the number
		/// of constraints. Welford's one-pass method is used. The number
		/// of constraints defaults to 1, applying Bessel's correction.
		/// A value of 0 may be used to compute the population standard deviation.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @param constraints The number of constraints, defaults to 1
		/// @return The standard deviation of the dataset
		template<typename Dataset>
		inline real stdev(const Dataset& data, unsigned int constraints = 1) {
			return sqrt(variance(data, constraints));
		}


		/// Compute the standard deviation of the mean given a dataset.
		/// Welford's one-pass method is used and Bessel's correction is applied.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The standard deviation of the mean
		template<typename Dataset>
		inline real stdom(const Dataset& X) {
			return sqrt(variance(X) / X.size());
		}


		/// Compute the relative error on a dataset using estimates of its mean
		/// and standard deviation, with the given number of constraints
		/// (defaults to 1 for Bessel's correction). The relative error is computed as
		/// \f$\epsilon_{rel} = \frac{\sigma}{\mu}\f$ and is not multiplied by 100.
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @param constraints The number of constraints for the estimators (defaults to 1)
		/// @return The standard relative error on the dataset
		template<typename Dataset>
		inline real standard_relative_error(const Dataset& X) {

			real x_mean = mean(X);

			if(abs(x_mean) < MACH_EPSILON) {
				TH_MATH_ERROR("standard_relative_error", x_mean, MathError::DivByZero);
				return nan();
			}

			return stdom(X) / abs(x_mean);
		}


		/// Compute the covariance between two datasets with the given number of constraints.
		/// The two datasets must have the same size.
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @param X The first dataset
		/// @param Y The second dataset
		/// @param constraints The number of constraints (defaults to 1 for Bessel's correction).
		/// @return The covariance between X and Y
		template<typename Dataset1, typename Dataset2>
		inline real covariance(
			const Dataset1& X, const Dataset2& Y, unsigned int constraints = 1) {

			if(X.size() != Y.size() || X.size() <= constraints) {
				TH_MATH_ERROR("covariance", X.size(), MathError::InvalidArgument);
				return nan();
			}

			real s = 0;
			real X_mean = mean(X);
			real Y_mean = mean(Y);

			for (unsigned int i = 0; i < X.size(); ++i)
				s += (X[i] - X_mean) * (Y[i] - Y_mean);

			return s / (X.size() - constraints);
		}


		/// Compute Pearson's correlation coefficient R between two datasets.
		/// The two datasets must have the same size.
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @param X The first dataset
		/// @param Y The second dataset
		/// @return The correlation coefficient computed using Pearson's formula
		template<typename Dataset1, typename Dataset2>
		inline real correlation_coefficient(
			const Dataset1& X, const Dataset2& Y) {

			return covariance(X, Y) / (stdev(X) * stdev(Y));
		}


		/// Compute the lag-n autocorrelation of a dataset as
		/// \f$\f$
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @param n The lag (defaults to lag-1)
		/// @return The lag-n autocorrelation of the given dataset
		template<typename Dataset>
		inline real autocorrelation(const Dataset& X, unsigned int n = 1) {

			if(X.size() < n) {
				TH_MATH_ERROR("autocorrelation", X.size(), MathError::InvalidArgument);
				return nan();
			}

			const real mu = mean(X);
			real num = 0;
			real den = square(X[0] - mu);

			for (unsigned int i = n; i < X.size(); ++i) {

				const real delta = X[i] - mu;
				num += delta * (X[i - n] - mu);
				den += delta * delta;
			}

			return num / den;
		}


		/// Compute the mean absolute deviation of a dataset as
		/// \f$\frac{\sum_{i = 1}^n |x_i - \hat \mu|}{n}\f$
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The mean absolute deviation of the dataset
		template<typename Dataset>
		inline real absolute_deviation(const Dataset& X) {

			real mu = mean(X);
			real res = 0;

			for (real x : X)
				res += abs(x - mu);

			return res / X.size();
		}


		/// Compute the skewness of a dataset as
		/// \f$\frac{\sum_{i=1}^n (\frac{x_i - \hat \mu}{\hat \sigma})^3}{n}\f$
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The skewness of the dataset
		template<typename Dataset>
		inline real skewness(const Dataset& X) {

			real mu, sigma;
			real res = 0;

			moments2(X, mu, sigma);
			sigma = sqrt(sigma);

			for (real x : X)
				res += cube((x - mu) / sigma);

			return res / X.size();
		}


		/// Compute the normalized kurtosis of a dataset as
		/// \f$\frac{\sum_{i=1}^n (\frac{x_i - \hat \mu}{\hat \sigma})^4}{n} - 3\f$
		///
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param X The dataset
		/// @return The normalized kurtosis of the dataset
		template<typename Dataset>
		inline real kurtosis(const Dataset& X) {

			real mu, sigma;
			real res = 0;

			moments2(X, mu, sigma);
			sigma = sqrt(sigma);

			for (real x : X)
				res += pow((x - mu) / sigma, 4);

			return (res / X.size()) - 3;
		}


		/// Compute the expectation value of a given function with respect
		/// to a Gaussian distribution with the given parameters.
		/// This function uses Gauss-Hermite quadrature to compute
		/// the integral \f$\int_{-\infty}^{+\infty} g(x) e^{-x^2} dx\f$
		///
		/// @tparam RealFunction A function or lambda representing a univariate real function
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
		/// @tparam Dataset Any type representing a dataset as a vector of values
		/// @param The data set to normalize
		/// @return The normalized data set
		template<typename Dataset>
		inline Dataset normalize_z_score(const Dataset& X) {

			real mu, sigma;
			moments2(X, mu, sigma);
			sigma = sqrt(sigma);

			return map([mu, sigma](real x) { return z_score(x, mu, sigma); }, X);
		}


		/// Compute the chi-square from the set of observed
		/// quantities, expected quantities and errors.
		/// The provided sets should all have the same size.
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @tparam Dataset3 Any type representing a dataset as a vector of values
		/// @param O The set of observed values
		/// @param E The set of expected values
		/// @param sigma The set of standard deviations on the observations
		/// @return The computed Chi-squared
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline real chi_square(
			const Dataset1& O, const Dataset2& E, const Dataset3& sigma) {

			if(O.size() != E.size() || E.size() != sigma.size()) {
				TH_MATH_ERROR("chi_square", E.size(), MathError::InvalidArgument);
				return nan();
			}

			real c_sqr = 0;

			for (unsigned int i = 0; i < O.size(); ++i) {

				if(abs(sigma[i]) < MACH_EPSILON) {
					TH_MATH_ERROR("chi_square", sigma[i], MathError::DivByZero);
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
				TH_MATH_ERROR("pvalue_chi_squared", ndf, MathError::InvalidArgument);
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
	  

		/// Compute the chi-square on a linear regression, as the sum of the squares
		/// of the residuals divided by the standard deviation.
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @tparam Dataset3 Any type representing a dataset as a vector of values
		///
		/// @param X A vector of the X values of the sample
		/// @param Y A vector of the Y values of the sample
		/// @param sigma The standard deviations of each point of the sample
		/// @param intercept The intercept of the linear model
		/// @param slope The slope of the linear model
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline real chi_square_linear(
			const Dataset1& X, const Dataset2& Y,
			const Dataset3& sigma, real intercept, real slope) {

			if(X.size() != Y.size() || X.size() != sigma.size()) {
				TH_MATH_ERROR(
					"chi_square_linear",
					X.size(), MathError::InvalidArgument);
				return nan();
			}

			real chi_squared = 0;
			for (unsigned int i = 0; i < X.size(); ++i) {

				if(abs(sigma[i]) <= MACH_EPSILON) {
					TH_MATH_ERROR("chi_square_linear", sigma[i], MathError::DivByZero);
					return nan();
				}

				chi_squared += square((Y[i] - intercept - slope * X[i]) / sigma[i]);
			}

			return chi_squared;
		}


		/// Compute the reduced chi-squared on a linear regression, computed as the usual
		/// chi-square (computed by chi_square_linear) divided by the number of degrees
		/// of freedom of the model (\f$N - 2\f$).
		///
		/// @tparam Dataset1 Any type representing a dataset as a vector of values
		/// @tparam Dataset2 Any type representing a dataset as a vector of values
		/// @tparam Dataset3 Any type representing a dataset as a vector of values
		///
		/// @param X A vector of the X values of the sample
		/// @param Y A vector of the Y values of the sample
		/// @param sigma The standard deviations of each point of the sample
		/// @param intercept The intercept of the linear model
		/// @param slope The slope of the linear model
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline real reduced_chi_square_linear(
			const Dataset1& X, const Dataset2& Y,
			const Dataset3& sigma, real intercept, real slope) {

			if(Y.size() <= 2) {
				TH_MATH_ERROR("reduced_chi_square_linear",
					Y.size(), MathError::InvalidArgument);
				return nan();
			}

			// Divide by degrees of freedom (N - 2)
			return chi_square_linear(X, Y, sigma, intercept, slope)
				/ (real) (Y.size() - 2);
		}
	}
}

#endif
