
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


namespace theoretica {


	/// @class histogram
	/// Histogram class with running statistics
	class histogram {
		public:

			/// Number of total data points
			size_t N;

			/// Bins
			std::vector<unsigned int> bins;

			/// Upper extreme of the interval to consider
			real range_max;

			/// Lower extreme of the interval to consider
			real range_min;
			
			/// Maximum value of the data
			real max_val;

			/// Minimum value of the data
			real min_val;

			/// Sum of the data points
			real sum;

			/// Square sum of the data points
			real sqr_sum;


			/// Default constructor
			histogram() :
			N(0), range_max(0), range_min(0),
			max_val(0), min_val(0), sum(0), sqr_sum(0) {}


			/// Construct from parameters
			histogram(unsigned int bin_count, real max, real min)
				: N(0), max_val(0), min_val(0), sum(0), sqr_sum(0) {

				bins.resize(bin_count);
				this->range_max = max;
				this->range_min = min;
			}


			/// Construct from data
			template<typename Dataset>
			histogram(const Dataset& data, unsigned int bin_count = 0) {

				range_max = max(data);
				range_min = min(data);
				max_val = range_max;
				min_val = range_min;
				sum = theoretica::sum(data);
				sqr_sum = sum_squares(data);
				N = data.size();

				// Default bin count is sqrt(N)
				bins.resize(
					bin_count ? bin_count : floor(sqrt(N))
				);

				// The histogram contains the data point
				// by construction
				for (size_t i = 0; i < N; ++i)
					bins[index(data[i])]++;
			}
		

			/// Insert a new data point inside the histogram.
			/// 
			/// @param x The value to insert
			inline void insert(real x) {

				if(x < range_min || x > range_max)
					return;

				bins[index(x)]++;
				N++;

				sum += x;
				sqr_sum += square(x);

				this->max_val = max(this->max_val, x);
				this->min_val = min(this->min_val, x);
			}


			/// Find the index corresponding to a given data point.
			///
			/// @param x The value to find the bin index of
			/// (must be between range_min and range_max)
			///
			/// @note This function does not check whether the
			/// value is between range_min and range_max.
			inline unsigned int index(real x) const {

				if(abs(x - range_max) < MACH_EPSILON)
					return bins.size() - 1;

				return floor(
					(x - range_min) / (range_max - range_min)
					* bins.size()
				);
			}


			// Statistics

			/// Get the number of points
			inline unsigned int number() {
				return N;
			}


			/// Get the mean value of the histogram data
			inline real mean() const {

				if(N == 0) {
					TH_MATH_ERROR("histogram::mean", N, DIV_BY_ZERO);
					return nan();
				}

				return sum / N;
			}


			/// Get the variance of the histogram data
			/// (with Bessel's correction)
			inline real variance() const {

				if(N == 0 || N == 1) {
					TH_MATH_ERROR("histogram::variance", N, DIV_BY_ZERO);
					return nan();
				}

				// TO-DO Implement Welford's method
				// to avoid catastrophic cancellation

				return (sqr_sum / N - square(sum / N)) * N / (real) (N - 1);
			}


			/// Get the standard deviation of the histogram data
			inline real stdev() const {
				return sqrt(variance());
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the histogram to string representation
			inline std::string to_string(
				const std::string& separator = " ") const {

				if(N == 0)
					return "";

				std::stringstream res;

				for (size_t i = 0; i < bins.size(); ++i) {

					res << (range_min + ((i + 1) / (real) bins.size())
							* (range_max - range_min))
						<< separator
						<< (bins[i] / (real) N)
						<< std::endl;
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

}

#endif
