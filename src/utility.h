
///
/// @file utility.h Utility functions to print or insert vec_buff data sets.
/// This header file is **not** automatically included by theoretica.h
///

#ifndef THEORETICA_UTILITY_H
#define THEORETICA_UTILITY_H

#include <iostream>
#include <string>
#include <algorithm>

#include "./core/vec_buff.h"
#include "./statistics/statistics.h"


namespace theoretica {


	/// Print a `vec_buff` data set to standard output
	void print_vec_buff(const vec_buff& v) {
		for(unsigned int i = 0; i < v.size(); i++) {
			std::cout << v[i] << std::endl;
		}
	}


	/// Print a `vec_buff` data set to standard output on a single row
	void print_vec_buff_row(const vec_buff& v) {
		std::cout << "{ ";
		for(unsigned int i = 0; i < v.size(); i++) {
			std::cout << v[i];

			if(i != v.size() - 1)
				std::cout << ", ";
		}
		std::cout << " }" << std::endl;
	}


	/// Insert a `vec_buff` data set from standard input
	void insert_data(vec_buff& data, std::string terminator) {

		std::string line;
		real value;

		while(true) {
			std::getline(std::cin, line);

			if(line == "")
				continue;

			if(line == terminator)
				break;

			std::replace(line.begin(), line.end(), ',', '.');

			try {
				value = std::stod(line);
			} catch(...) {
				std::cout << "Input conversion error" << std::endl;
				value = nan();
			}

			if(!is_nan(value))
				data.emplace_back(value);
		}
	}


	/// Print common statistical information about a `vec_buff` data set
	void print_sample_stats(const vec_buff& X) {
		std::cout << "N = " << X.size() << std::endl;
		std::cout << "Mean: " << mean(X) << std::endl;
		std::cout << "Variance: " << sample_variance(X) << std::endl;
		std::cout << "Standard Deviation: " << smpl_stdev(X) << std::endl;
		std::cout << "Relative Error: " << sample_standard_relative_error(X) * 100 << "%" << std::endl;
		std::cout << "Mean Standard Deviation: " << smpl_stdom(X) << std::endl;
		std::cout << "Chi-Square (Sigma): " << chi_square_sigma(X) << std::endl;
	}

}

#endif
