#ifndef UROBORO_STATISTICS_H
#define UROBORO_STATISTICS_H

#include "./constants.h"
#include "./vec.h"
#include "./vec_buff.h"
#include "./common.h"

namespace uroboro {

	// Calculate the mean of a set of values
	template<unsigned int N>
	inline real mean(const vec<N>& dataset) {
		real sum = 0;
		for (int i = 0; i < N; ++i) {
			sum += dataset.data[i];
		}
		return sum / (real) N;
	}

	inline real mean(const vec_buff& dataset) {
		real sum = 0;
		for (auto i : dataset) {
			sum += i;
		}
		return sum / (real) dataset.size();
	}


	// Calculate the weighted mean of a set of values
	// <dataset> and <weights> must have the same size
	template<unsigned int N>
	inline real weighted_mean(const vec<N>& dataset, const vec<N>& weights) {
		real sum = 0;
		real weight_sum = 0;
		for (int i = 0; i < N; ++i) {
			sum += dataset.data[i] * weights.data[i];
			weight_sum += weights.data[i];
		}
		return sum / weight_sum;
	}

	inline real weighted_mean(const vec_buff& dataset, const vec_buff& weights) {

		if(dataset.size() != weights.size())
			// throw ...
			return 0;

		real sum = 0;
		real weight_sum = 0;
		for (int i = 0; i < dataset.size(); ++i) {
			sum += dataset[i] * weights[i];
			weight_sum += weights[i];
		}
		return sum / weight_sum;
	}


	// Calculate the variance of a population
	template<unsigned int N>
	inline real variance(const vec<N>& dataset) {
		real sum = 0;
		real Xm = mean(dataset);
		real diff = 0;
		for (int i = 0; i < N; ++i) {
			diff = dataset.data[i] - Xm;
			sum += diff * diff;
		}
		return sum / (real) N;
	}

	inline real variance(const vec_buff& dataset) {
		real sum = 0;
		real Xm = mean(dataset);
		real diff = 0;
		for (auto i : dataset) {
			diff = i - Xm;
			sum += diff * diff;
		}
		return sum / (real) dataset.size();
	}


	// Calculate the variance of a sample
	// This function uses Bessel correction
	template<unsigned int N>
	inline real sample_variance(const vec<N>& dataset) {
		real sum = 0;
		real Xm = mean(dataset);
		real diff = 0;
		for (int i = 0; i < N; ++i) {
			diff = dataset.data[i] - Xm;
			sum += diff * diff;
		}
		return sum / ((real) N - 1);
	}

	inline real sample_variance(const vec_buff& dataset) {
		real sum = 0;
		real Xm = mean(dataset);
		real diff = 0;
		for (auto i : dataset) {
			diff = i - Xm;
			sum += diff * diff;
		}
		return sum / (real) (dataset.size() - 1);
	}


	// Calculate the standard deviation of a population
	template<unsigned int N>
	inline real standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(variance(dataset));
	}

	inline real standard_deviation(const vec_buff& dataset) {
		return uroboro::sqrt(variance(dataset));
	}


	// Calculate the standard deviation of a sample
	template<unsigned int N>
	inline real sample_standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(sample_variance(dataset));
	}

	inline real sample_standard_deviation(const vec_buff& dataset) {
		return uroboro::sqrt(sample_variance(dataset));
	}


	// Calculate the standard deviation on the mean of a set of measures
	// Bessel correction is used in the calculation of variance
	template<unsigned int N>
	inline real mean_standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(sample_variance(dataset)) / uroboro::sqrt(dataset.size());
	}

	inline real mean_standard_deviation(const vec_buff& dataset) {
		return uroboro::sqrt(sample_variance(dataset)) / uroboro::sqrt(dataset.size());
	}


	// Calculate the covariance of two sets of measures
	template<unsigned int N>
	inline real covariance(const vec<N>& X, const vec<N>& Y) {
		real sum = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);
		for (int i = 0; i < N; ++i) {
			sum += (X.data[i] - X_mean) * (Y.data[i] - Y_mean);
		}
		return sum / (real) N;
	}

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
	template<unsigned int N>
	inline real sample_covariance(const vec<N>& X, const vec<N>& Y) {
		real sum = 0;
		real X_mean = mean(X);
		real Y_mean = mean(Y);
		for (int i = 0; i < N; ++i) {
			sum += (X.data[i] - X_mean) * (Y.data[i] - Y_mean);
		}
		return sum / (real) (N - 1);
	}

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
		return sum / (real) (X.size() - 1);
	}

}

#endif
