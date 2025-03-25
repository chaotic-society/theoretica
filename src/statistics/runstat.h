
///
/// @file runstat.h Classes for running statistics computation.
///

#ifndef THEORETICA_RUNSTAT_H
#define THEORETICA_RUNSTAT_H

#include "../core/constants.h"


namespace theoretica {
namespace stats {


	// A running statistics object is defined as a class which
	// provides, at least, a method insert() to add values to the sample,
	// a method get() to get the computed running statistic and
	// a method clear() to remove all points from the sample.
	// The class may also provide specific methods to obtain other estimates.


	///
	/// @class runstat_moments2_t
	/// A running statistics class which computes the mean and variance of
	/// the provided data points. By default, the Type of variables is real
	/// and the runstat_moments2 typedef is available, but it may be also
	/// used for vectors or any type which has arithmetic operators.
	/// 
	template <
		typename Type = real
	>
	class runstat_moments2_t {
	private:

		/// Running average
		Type average {0.0};


		/// Running total sum of squares
		Type tss {0.0};


		/// Size of the sample
		unsigned int sample_size {0};


	public:
		
		/// Default constructor
		runstat_moments2_t() {}


		/// Insert a new data value in the running statistics sample computation.
		///
		/// @param x The value to insert.
		inline runstat_moments2_t<Type>& insert(Type x) {

			// Update statistics using Welford's method
			const Type tmp = average;
			average = tmp + (x - tmp) / (sample_size + 1);
			tss = tss + (x - tmp) * (x - average);
			sample_size++;

			return *this;
		}


		/// Get the estimated mean of the sample.
		inline Type mean() {
			return average;
		}


		/// Get the variance of the sample.
		inline Type variance() {

			if (sample_size <= 1) {
				TH_MATH_ERROR(
					"stats::runstat_moments2_t::variance", sample_size, INVALID_ARGUMENT);
				return Type(nan());
			}

			return tss * (1.0 / (sample_size - 1));
		}


		/// Get the number of values in the sample.
		inline unsigned int number() {
			return sample_size;
		}


		/// Get the estimated statistic, as a vector of mean and variance.
		inline Type get() {
			return variance();
		}


		/// Clear the stored statistics and internal counters.
		inline runstat_moments2_t<Type>& clear() {
			
			average = 0.0;
			tss = 0.0;
			sample_size = 0;

			return *this;
		}
	};


	/// runstat_moments2 type for real random variables.
	using runstat_moments2 = runstat_moments2_t<real>;
	

}}

#endif
