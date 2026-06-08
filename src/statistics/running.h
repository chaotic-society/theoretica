
///
/// @file running.h Classes for running statistics.
///

#ifndef THEORETICA_RUNNING_STATISTICS_H
#define THEORETICA_RUNNING_STATISTICS_H

#include "../core/constants.h"
#include "../core/core_traits.h"


namespace theoretica {
namespace stats {


	/// @class RunningMoments2
	/// A running statistics class which computes the mean and variance of
	/// the provided data points. By default, the Type of variables is real,
	/// but it may be also used for vectors or any type which has arithmetic operators.
	template <typename Type = real>
	class RunningMoments2 {
	private:

		/// Running average
		Type average {};

		/// Running total sum of squares
		Type tss {};

		/// Size of the sample
		size_t sample_size {0};


	public:
		
		/// Default constructor
		RunningMoments2() {}


		/// Insert a new data value in the running statistics sample computation.
		/// Non-vector overload (scalars or non-vector types).
		template <
			typename ScalarType = Type,
			disable_vector<ScalarType> = true
		>
		inline RunningMoments2<Type>& insert(const ScalarType& x) {

			// Update statistics using Welford's method
			const Type tmp = average;
			average = tmp + (x - tmp) / (sample_size + 1);
			tss = tss + (x - tmp) * (x - average);
			sample_size++;

			return *this;
		}


		/// Insert a new data value in the running statistics sample computation.
		/// Vector overload: initialize on first insert then update per-component.
		template <
			typename VectorType = Type,
			enable_vector<VectorType> = true
		>
		inline RunningMoments2<Type>& insert(const VectorType& x) {

			// Initialize average and tss appropriately
			if(sample_size == 0) {

				average = x;
				tss.resize(x.size());

				// Initialize tss to zero per-component
				for (unsigned int i = 0; i < tss.size(); ++i)
					tss[i] = 0.0;

				sample_size = 1;
				return *this;
			}

			const Type tmp = average;
			for (unsigned int i = 0; i < x.size(); ++i) {
				average[i] = tmp[i] + (x[i] - tmp[i]) / (sample_size + 1);
				tss[i] = tss[i] + (x[i] - tmp[i]) * (x[i] - average[i]);
			}

			sample_size++;

			return *this;
		}


		/// Get the estimated mean of the sample.
		inline Type mean() const {
			return average;
		}


		/// Get the variance of the sample.
		inline Type variance() const {

			if (sample_size <= 1) {
				TH_MATH_ERROR("stats::RunningMoments2::variance", sample_size, MathError::ImpossibleOperation);
				return make_error<Type>();
			}

			return tss * (1.0 / (sample_size - 1));
		}


		/// Get the number of values in the sample.
		inline size_t number() const {
			return sample_size;
		}


		/// Get the estimated statistic, as a 2D vector of mean and variance.
		inline std::pair<Type, Type> get() const {
			return { mean(), variance() };
		}


		/// Clear the stored statistics and internal counters.
		inline RunningMoments2<Type>& clear() {
            
			average = Type{};
			tss = Type{};
			sample_size = 0;

			return *this;
		}


		/// Combine the mean and variance computed
		/// by another RunningMoments2 object with the current one,
		/// modifying its content. Chan's parallel algorithm is used for the combination.
		///
		/// @param other The other RunningMoments2 object to combine with the current one.
		/// @return A reference to the current object, with the combined statistics.
		inline RunningMoments2<Type>& combine(const RunningMoments2<Type>& other) {

			if(other.sample_size == 0)
				return *this;

			if(sample_size == 0) {
				*this = other;
				return *this;
			}

			const Type delta = other.mean() - mean();
			const size_t total_size = sample_size + other.sample_size;
			const real ratio = real(other.sample_size) / total_size;

			average = mean() + delta * ratio;
			tss = tss + other.tss + delta * delta * (real(sample_size) * ratio);
			sample_size = total_size;

			return *this;
		}
	};

}}

#endif