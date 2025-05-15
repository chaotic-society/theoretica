
///
/// @file histogram.h Histogram class
///

#ifndef THEORETICA_HISTOGRAM_H
#define THEORETICA_HISTOGRAM_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include <vector>
#include "../core/real_analysis.h"
#include "../core/dataset.h"
#include "./statistics.h"


namespace theoretica {


	/// @class histogram
	/// Histogram class with running statistics, can be constructed
	/// from the parameters of the bins or from a dataset. Elements
	/// are inserted one by one, updating the running statistics
	/// for the TSS, mean, maximum and minimum on each step.
	class histogram {
		private:

			/// Number of total data points
			size_t N {0};

			/// Bins
			std::vector<unsigned int> bin_counts;

			/// Upper extreme of the interval to consider
			real range_max;

			/// Lower extreme of the interval to consider
			real range_min;
			
			/// Maximum value of the data
			real value_max;

			/// Minimum value of the data
			real value_min;

			/// Running average
			real run_average;

			/// Running total sum of squares
			real run_tss;

		public:

			/// Construct the histogram from the number of bins and the range.
			///
			/// The histogram is initialized from the arguments, without specifying
			/// any data points, which need to be added with insert().
			///
			/// @param bin_count The number of bins
			/// @param range_min The lower bound of the range
			/// @param range_max The upper bound of the range
			histogram(unsigned int bin_count, real range_min, real range_max)
				: N(0), value_max(-inf()), value_min(-inf()), run_average(0), run_tss(0) {

				bin_counts.resize(bin_count);
				this->range_max = range_max;
				this->range_min = range_min;
			}


			/// Construct the histogram from a set of data points, with the given
			/// number of bins. If the number of bins is not specified,
			/// it defaults to \f$[\sqrt{N}]\f$.
			///
			/// @param data The set of data points
			/// @param bin_count The number of bins
			/// (defaults to the square root of the number of points)
			template<typename Dataset, enable_vector<Dataset> = true>
			histogram(const Dataset& data, unsigned int bin_count = 0) {

				range_max = theoretica::max(data);
				range_min = theoretica::min(data);
				value_max = range_max;
				value_min = range_min;
				N = data.size();

				// Compute mean and TSS
				run_average = stats::mean(data);
				run_tss = stats::total_sum_squares(data);

				// Default bin count is sqrt(N)
				bin_counts.resize(
					bin_count ? bin_count : floor(sqrt(N))
				);

				// The histogram contains all the data points by construction
				for (size_t i = 0; i < N; ++i)
					bin_counts[index(data[i])]++;
			}


			/// Insert a new data point inside the histogram, updating
			/// the running statistics and the corresponding bin.
			///
			/// @param x The value to insert
			inline void insert(real x) {

				if(x < range_min || x > range_max)
					return;

				// Update average and TSS using Welford's method
				const real tmp = run_average;
				run_average = tmp + (x - tmp) / (N + 1);
				run_tss += (x - tmp) * (x - run_average);

				value_max = value_max < x ? x : value_max;
				value_min = value_min > x ? x : value_min;

				bin_counts[index(x)]++;
				N++;
			}


			/// Find the bin index corresponding to a given data point.
			///
			/// @param x The value to find the bin index of
			/// (must be between range_min and range_max)
			///
			/// @note This function does not check whether the
			/// value is between range_min and range_max,
			/// so care should be taken to use only valid inputs.
			inline unsigned int index(real x) const {

				if(abs(x - range_max) < MACH_EPSILON)
					return bin_counts.size() - 1;

				return floor(
					(x - range_min) / (range_max - range_min)
					* bin_counts.size()
				);
			}


			// Statistical functions


			/// Get the number of data points inside the histogram
			///
			/// @return The number of data points which have been
			/// added to the histogram.
			inline unsigned int number() const {
				return N;
			}


			/// Get a vector containing the bin counts of each bin.
			///
			/// @note The bins cannot be directly modified,
			/// new elements must be added using insert().
			///
			/// @return A vector containing the number of elements in each bin.
			inline std::vector<unsigned int> bins() const {
				return bin_counts;
			}


			/// Get the biggest data point of the histogram.
			///
			/// @return The maximum value of all elements.
			inline real max() const {
				return value_max;
			}


			/// Get the smallest data point of the histogram.
			///
			/// @return The minimum value of all elements.
			inline real min() const {
				return value_min;
			}


			/// Get the mean value of the histogram data.
			///
			/// @return The running mean of all elements of the histogram.
			inline real mean() const {
				return run_average;
			}


			/// Get the total sum of squares (TSS) computed
			/// using Welford's one-pass method.
			///
			/// @return The total sum of squares of all elements of the histogram.
			inline real tss() const {
				return run_tss;
			}


			// Operators


			/// Evaluate the histogram like a step function
			/// which is zero outside the range of the histogram.
			///
			/// @param x The point to evaluate the histogram function at
			/// @return The value of the histogram function at x
			inline real operator()(real x) {

				if (x < range_min || x > range_max)
					return 0.0;

				return bin_counts[index(x)];
			}


			/// Get the number of elements in the i-th bin.
			///
			/// @param i The index of the bin
			/// @return The number of elements in the i-th bin
			inline unsigned int operator[](unsigned int i) const {
				return bin_counts[i];
			}


			// TO-DO Cumulative Distribution Function


#ifndef THEORETICA_NO_PRINT

			/// Convert the histogram to string representation
			///
			/// @param separator The string to print between row elements
			/// @param normalized Whether to normalize the bin counts as a frequency
			/// (defaults to true).
			/// @param lower_extreme Whether to print the lower extreme of the bins
			/// or use the mid point (defaults to false, using mid points).
			/// @return A string representing the histogram, ready to plot.
			inline std::string to_string(
				const std::string& separator = " ",
				bool normalized = true,
				bool lower_extreme = false) const {

				if(N == 0)
					return "";

				std::stringstream res;
				const real width = abs(range_max - range_min) / bin_counts.size();
				real mult = 0.5;

				if (lower_extreme)
					mult = 0.0;

				for (size_t i = 0; i < bin_counts.size(); ++i) {

					res << (range_min + (i + mult) * width) << separator;

					if (normalized)
						res << (bin_counts[i] / (real) N) << std::endl;
					else
						res << bin_counts[i] << std::endl;
				}

				return res.str();
			}


			/// Convert the histogram to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the histogram in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const histogram& obj) {
				return out << obj.to_string();
			}

#endif

	};


	// Statistical functions over elements of a histogram


	namespace stats {


		/// Compute the mean of the values of a histogram.
		inline real mean(const histogram& h) {
			return h.mean();
		}


		/// Compute the total sum of squares of the values of the histogram.
		inline real tss(const histogram& h) {
			return h.tss();
		}


		/// Compute the variance of the values of a histogram.
		inline real variance(const histogram& h) {

			if (h.number() <= 1) {
				TH_MATH_ERROR("variance", h.number(), DIV_BY_ZERO);
				return nan();
			}

			return h.tss() / (h.number() - 1);
		}


		/// Compute the standard deviation of the values of a histogram.
		inline real stdev(const histogram& h) {
			return sqrt(variance(h));
		}
	}


	/// Compute the maximum value of the elements of a histogram.
	inline real max(const histogram& h) {
		return h.max();
	}


	/// Compute the minimum value of the elements of a histogram.
	inline real min(const histogram& h) {
		return h.min();
	}

}

#endif
