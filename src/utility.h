#ifndef UROBORO_UTILITY_H
#define UROBORO_UTILITY_H

#include <iostream>
#include <string>
#include <algorithm>

#include "./vec_buff.h"
#include "./statistics/statistics.h"


namespace uroboro {


	void print_vec_buff(const vec_buff& v) {
		for(int i = 0; i < v.size(); i++) {
			std::cout << v[i] << std::endl;
		}
	}


	void print_vec_buff_row(const vec_buff& v) {
		std::cout << "{ ";
		for(int i = 0; i < v.size(); i++) {
			std::cout << v[i];

			if(i != v.size() - 1)
				std::cout << ", ";
		}
		std::cout << " }" << std::endl;
	}


	void insert_data(vec_buff& data, std::string terminator) {

        std::string line;
        real value;

        while(true) {
            std::getline(std::cin, line);

            if(line == terminator)
                break;

            std::replace(line.begin(), line.end(), ',', '.');

            try {
                value = std::stod(line);
            } catch(...) {
                std::cout << "Input conversion error" << std::endl;
            }

            data.emplace_back(value);
        }
	}


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
