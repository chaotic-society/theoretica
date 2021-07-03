#ifndef UROBORO_STATISTICS_H
#define UROBORO_STATISTICS_H

#include "./constants.h"
#include "./vec.h"
#include "./vec_buff.h"
#include "./common.h"


namespace uroboro {


	// Calculate the mean of a set of values
	inline real mean(const vec_buff& data) {

		// Sum of x_i / N
		return sum(data) / (real) data.size();
	}


	// Calculate the weighted mean of a set of values
	// <data> and <weights> must have the same size
	inline real weighted_mean(const vec_buff& data, const vec_buff& weights) {

		// Sum of x_i * w_i / Sum of w_i
		return product_sum(data, weights) / sum(weights);
	}


	// Propagate error of a sum of values
	// as sqrt(sigma_x^2 + sigma_y^2 + ...)
	inline real propagate_sum(const vec_buff& sigma) {

		return uroboro::sqrt(sum_squares(sigma));
	}


	// Propagate error of a product (or quotient) of values
	// as sqrt((sigma_x / x_mean)^2 + (sigma_y / y_mean)^2 + ...)
	// The result is the propagated relative error
	inline real propagate_product(const vec_buff& sigma, const vec_buff& mean) {

		if(sigma.size() != mean.size())
			// throw...
			return 0;

		// Calculate sum of squares of (i_sigma / i_mean)
		real sum = 0;
		for (int i = 0; i < sigma.size(); ++i) {
			sum += square(sigma[i] / uroboro::abs(mean[i]));
		}

		return uroboro::sqrt(sum);
	}


	// Total sum of squares (TSS)
	// Calculated as sum(square(x_i - x_mean))
	inline real total_sum_squares(const vec_buff& X) {

		if(!X.size())
			//throw...
			return 0;

		real tss = 0;
		real x_m = mean(X);
		for (auto x : X) {
			tss += square(x - x_m);
		}

		return tss;
	}

	// Total sum of squares (TSS)
	// Calculated as sum(square(x_i - x_mean))
	inline real tss(const vec_buff& X) {
		return total_sum_squares(X);
	}


	// Calculate the variance of a population
	inline real variance(const vec_buff& data) {

		if(!data.size())
			//throw...
			return 0;

		return total_sum_squares(data) / (real) data.size();
	}


	// Calculate the variance of a sample
	// This function uses Bessel correction
	inline real sample_variance(const vec_buff& data) {

		if(!data.size())
			//throw...
			return 0;

		// Bessel correction (N - 1)
		return total_sum_squares(data) / (real) (data.size() - 1);
	}


	// Calculate the standard deviation of a population
	inline real standard_deviation(const vec_buff& data) {
		return uroboro::sqrt(variance(data));
	}


	// Calculate the standard deviation of a population
	inline real stdev(const vec_buff& data) {
		return standard_deviation(data);
	}


	// Calculate the standard deviation of a sample
	inline real sample_standard_deviation(const vec_buff& data) {
		return uroboro::sqrt(sample_variance(data));
	}


	// Calculate the standard deviation of a sample
	inline real smpl_stdev(const vec_buff& data) {
		return sample_standard_deviation(data);
	}


	// Calculate the relative error on a population measure
	// using standard deviation
	inline real standard_relative_error(const vec_buff& X) {
		return standard_deviation(X) / uroboro::abs(mean(X));
	}


	// Calculate the relative error on a sample measure
	// using standard deviation
	inline real sample_standard_relative_error(const vec_buff& X) {
		return sample_standard_deviation(X) / uroboro::abs(mean(X));
	}


	// Calculate the standard deviation on the mean of a set of values
	inline real mean_standard_deviation(const vec_buff& data) {
		return uroboro::sqrt(variance(data)) / uroboro::sqrt(data.size());
	}


	// Calculate the standard deviation on the mean of a set of values
	inline real stdom(const vec_buff& data) {
		return mean_standard_deviation(data);
	}


	// Calculate the standard deviation on the mean of a set of measures
	// Bessel correction is used in the calculation of variance
	inline real sample_mean_standard_deviation(const vec_buff& data) {
		return uroboro::sqrt(sample_variance(data)) / uroboro::sqrt(data.size());
	}


	// Calculate the standard deviation on the mean of a set of measures
	// Bessel correction is used in the calculation of variance
	inline real smpl_stdom(const vec_buff& data) {
		return sample_mean_standard_deviation(data);
	}


	// Calculate the covariance of two sets of measures
	inline real covariance(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw ...
			return 0;

		real sum = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);
		for (int i = 0; i < X.size(); ++i) {
			sum += (X[i] - X_mean) * (Y[i] - Y_mean);
		}
		return sum / (real) X.size();
	}


	// Calculate the covariance between two sets of sample measures
	// This function uses Bessel correction
	inline real sample_covariance(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw ...
			return 0;

		real sum = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);
		for (int i = 0; i < X.size(); ++i) {
			sum += (X[i] - X_mean) * (Y[i] - Y_mean);
		}

		// Bessel correction (N - 1)
		return sum / (real) (X.size() - 1);
	}


	// Pearson's correlation coefficient R for a population
	inline real correlation_coefficient(const vec_buff& X, const vec_buff& Y) {
		return covariance(X, Y) / (stdev(X) * stdev(Y));
	}


	// Pearson's correlation coefficient r for a sample
	inline real sample_correlation_coefficient(const vec_buff& X, const vec_buff& Y) {
		return sample_covariance(X, Y) / (smpl_stdev(X) * smpl_stdev(Y));
	}


	// Gaussian Distribution function
	inline real gaussian_distribution(real x, real X, real sigma) {

		return (1.0 / (sigma *
			uroboro::sqrt(2 * PI))) * uroboro::exp(-square(x - X) / (2 * square(sigma)));
	}


	// Gaussian Distribution function calculated on a sample of measures
	inline real gaussian_distribution(real x, const vec_buff& data) {

		return gaussian_distribution(
			x, mean(data),
			sample_mean_standard_deviation(data));
	}


	inline real binomial_distribution(unsigned int nu, unsigned int n, real p) {
		return binomial_coeff(n, nu) *
			uroboro::pow(p, nu) * uroboro::pow(1 - p, n - nu);
	}


	// TO-DO
	// gaussian distribution probability inside (t * sigma)
	// erf


	// Normal distribution chi-square with 4 intervals
	// calculated on a sample of measures
	inline real chi_square_sigma(const vec_buff& X) {

		if(!X.size())
			//throw...
			return 0;

		const unsigned int N = X.size();

		unsigned int Ok_1 = 0;
		unsigned int Ok_2 = 0;
		unsigned int Ok_3 = 0;
		unsigned int Ok_4 = 0;

		real M = mean(X);
		real sigma = sample_standard_deviation(X);

		for (auto x : X) {
			if(x < M - sigma)
				Ok_1++;
			else if((x > M - sigma) && (x < M))
				Ok_2++;
			else if((x > M) && (x < M + sigma))
				Ok_3++;
			else if(x > M + sigma)
				Ok_4++;
		}

		// Sum of (Ok - Ek)^2 / Ek
		// where Ek = N * Pk
		real chi_2 = ((Ok_1 - (N * 0.16)) * (Ok_1 - (N * 0.16))) / (N * 0.16);
		chi_2 += ((Ok_2 - (N * 0.34)) * (Ok_2 - (N * 0.34))) / (N * 0.34);
		chi_2 += ((Ok_3 - (N * 0.34)) * (Ok_3 - (N * 0.34))) / (N * 0.34);
		chi_2 += ((Ok_4 - (N * 0.16)) * (Ok_4 - (N * 0.16))) / (N * 0.16);

		return chi_2;
	}


	// Calculate the intercept of the minimum squares linearization of X and Y
	inline real least_squares_linear_intercept(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw...
			return 0;

		real Delta = X.size() * product_sum(X, X) - square(sum(X));
		real A = (product_sum(X, X) * sum(Y) - sum(X) * product_sum(X, Y)) / Delta;

		return A;
	}


	// Calculate the intercept of the minimum squares linearization of X and Y
	inline real lst_sqrs_lin_intercept(const vec_buff& X, const vec_buff& Y) {
		return least_squares_linear_intercept(X, Y);
	}


	// Calculate the slope of the minimum squares linearization of X and Y
	inline real least_squares_linear_slope(const vec_buff& X, const vec_buff& Y) {

		if(X.size() != Y.size())
			// throw...
			return 0;

		real Delta = X.size() * product_sum(X, X) - square(sum(X));
		real B = (X.size() * product_sum(X, Y) - sum(X) * sum(Y)) / Delta;

		return B;
	}


	// Calculate the slope of the minimum squares linearization of X and Y
	inline real lst_sqrs_lin_slope(const vec_buff& X, const vec_buff& Y) {
		return least_squares_linear_slope(X, Y);
	}


	// Calculate the error of the minimum squares linearization of a sample
	real least_squares_linear_error(const vec_buff& X, const vec_buff& Y,
		real intercept, real slope) {

		if(X.size() != Y.size())
			// throw...
			return 0;

		real err = 0;
		for (int i = 0; i < X.size(); ++i) {
			err += square(Y[i] - intercept - slope * X[i]);
		}

		// Correction by degrees of freedom (N - 2)
		return uroboro::sqrt(err / (real) (X.size() - 2));
	}


	// Calculate the error of the minimum squares linearization of a sample
	inline real lst_sqrs_lin_error(const vec_buff& X, const vec_buff& Y,
		real intercept, real slope) {
		return least_squares_linear_error(X, Y, intercept, slope);
	}


	// Calculate the chi-square on a linearization
	real chi_square_linearization(const vec_buff& X, const vec_buff& Y, const vec_buff& sigma,
		real intercept, real slope) {

		if(X.size() != Y.size() || X.size() != sigma.size())
			// throw...
			return 0;

		real chi_squared = 0;
		for (int i = 0; i < X.size(); ++i) {
			chi_squared += square((Y[i] - intercept - slope * X[i]) / sigma[i]);
		}

		return chi_squared;
	}


	// Calculate the reduced chi-squared on a linearization
	real reduced_chi_square_linearization(const vec_buff& X, const vec_buff& Y, const vec_buff& sigma,
		real intercept, real slope) {

		// Divide by degrees of freedom (N - 2)
		return chi_square_linearization(X, Y, sigma, intercept, slope)
			/ (real) (Y.size() - 2);
	}


	// Calculate the intercept of the weighted minimum squares linearization of X and Y
	inline real least_squares_weighted_linear_intercept(const vec_buff& X,
		const vec_buff& Y, const vec_buff& W) {

		if(X.size() != Y.size() || X.size() != W.size())
			// throw...
			return 0;

		real Delta = sum(W) * product_sum(X, X, W) - square(product_sum(X, W));

		real A = (product_sum(X, X, W) * product_sum(Y, W) -
			product_sum(X, W) * product_sum(X, Y, W)) / Delta;

		return A;
	}


	// Calculate the intercept of the weighted minimum squares linearization of X and Y
	inline real lst_sqrs_weight_lin_intercept(const vec_buff& X, const vec_buff& Y, const vec_buff& W) {
		return least_squares_weighted_linear_intercept(X, Y, W);
	}


	// Calculate the slope of the weighted minimum squares linearization of X and Y
	inline real least_squares_weighted_linear_slope(const vec_buff& X,
		const vec_buff& Y, const vec_buff& W) {

		if(X.size() != Y.size() || X.size() != W.size())
			// throw...
			return 0;

		real Delta = sum(W) * product_sum(X, X, W) - square(product_sum(X, W));

		real B = (sum(W) * product_sum(X, Y, W) -
			product_sum(X, W) * product_sum(Y, W)) / Delta;

		return B;
	}


	// Calculate the slope of the weighted minimum squares linearization of X and Y
	inline real lst_sqrs_weight_lin_slope(const vec_buff& X, const vec_buff& Y, const vec_buff& W) {
		return least_squares_weighted_linear_slope(X, Y, W);
	}


}

#endif
