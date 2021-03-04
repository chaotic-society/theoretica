#ifndef UROBORO_STATISTICS_H
#define UROBORO_STATISTICS_H

#include "./vec.h"
#include "./common.h"

namespace uroboro {

	template<unsigned int N>
	inline real mean(const vec<N>& dataset) {
		real sum = 0;
		for (int i = 0; i < N; ++i) {
			sum += dataset.data[i];
		}
		return sum / (real) N;
	}

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

	template<unsigned int N>
	inline real standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(variance(dataset));
	}

	template<unsigned int N>
	inline real sample_standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(sample_variance(dataset));
	}

	template<unsigned int N>
	inline real mean_standard_deviation(const vec<N>& dataset) {
		return uroboro::sqrt(sample_variance(dataset)) / uroboro::sqrt(N);
	}

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

}

#endif
