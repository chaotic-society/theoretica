
///
/// @file bootstrap.h Statistical bootstrap methods.
///

#ifndef THEORETICA_BOOTSTRAP_H
#define THEORETICA_BOOTSTRAP_H

#include "./runstat.h"


namespace theoretica {
namespace stats {


	/// Monte Carlo Boostrap method to compute an arbitrary statistic from a sample.
	///
	/// @tparam RunStat The type of running statistics class
	/// to use (e.g. stats::runstat_moments).
	/// @tparam RunStatResult The type of running statistics class used on
	/// the resample statistics. By default, the mean and variance of the
	/// statistic estimated over the resamples are computed. The results are then
	/// returned using the class' get() method.
	/// @tparam Type The type of the statistic estimated (defaults to real).
	/// @tparam Dataset The type used to store the sample dataset.
	/// @param x The sample dataset.
	/// @param g An already initialized pseudorandom number generator.
	/// @param n The number of samples to generate.
	/// @return An object containing the mean and variance of the estimated 
	///
	/// The RunStat type is expected to be a running statistics class which takes in
	/// values of the type of the elements of Dataset with an insert() method
	/// and estimates statistics which are returned as a real value or vector using get().
	template <
		typename RunStat,
		typename Type = real,
		typename RunStatResult = runstat_moments2_t<Type>,
		typename Dataset
	>
	inline auto bootstrap(Dataset x, PRNG& g, unsigned int n = 1000) {

		// Compute the desired statistics of the estimator
		RunStatResult resample_stats {};

		if (!x.size()) {
			TH_MATH_ERROR("stats::bootstrap", x.size(), INVALID_ARGUMENT);
			return resample_stats;
		}
		
		// Generate n samples
		for (int i = 0; i < n; ++i) {

			RunStat stat = RunStat();

			// Generate samples with the same size as the original
			for (int j = 0; j < x.size(); ++j) {

				// Select a random value from the sample
				stat.add(x[g() % x.size()]);
			}

			resample_stats.insert(stat.get());
		}

		return resample_stats;
	}


	/// Monte Carlo Boostrap method to compute an arbitrary statistic from a sample.
	///
	/// @tparam RunStatResult The type of running statistics class used on
	/// the resample statistics. By default, the mean and variance of the
	/// statistic estimated over the resamples are computed. The results are then
	/// returned using the class' get() method.
	/// @tparam Dataset The type used to store the sample dataset.
	/// @tparam Estimator The type of the function used to estimate the statistic.
	/// @param x The sample dataset.
	/// @param estimate A function which computes the statistic of interest.
	/// @param g An already initialized pseudorandom number generator.
	/// @param n The number of samples to generate (defaults to 1000).
	/// @return An object containing the mean and variance of the estimated 
	///
	/// The RunStatResult type is expected to be a running statistics class which takes in
	/// values of the type of the elements of Dataset with an insert() method
	/// and estimates statistics which are returned as a real value or vector using get().
	/// The Dataset type is expected to provide a resize() method to change the size of
	/// the contained (such as std::vector or th::vec).
	template <
		typename RunStatResult = runstat_moments2,
		typename Dataset,
		typename Estimator = std::function<real(Dataset)>
	>
	auto bootstrap(Dataset x, Estimator estimate, PRNG& g, unsigned int n = 1000) {

		// Resampled vector
		Dataset resample;
		resample.resize(x.size());

		RunStatResult resample_stats {};

		// Construct n samples
		for (size_t i = 0; i < n; ++i) {

			for (size_t j = 0; j < x.size(); ++j)
				resample[j] = x[g() % x.size()];

			resample_stats.insert(estimate(resample));
		}

		// Compute the variance of the estimator
		return resample_stats.get();
	}

}}

#endif
