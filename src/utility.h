#ifndef UROBORO_UTILITY_H
#define UROBORO_UTILITY_H
#include <iostream>
#include <string>
#include <algorithm>

#include "./algebra/vec.h"
#include "./vec_buff.h"
#include "./algebra/mat.h"
#include "./complex/complex.h"
#include "./complex/quat.h"
#include "./complex/phasor.h"
#include "./statistics/statistics.h"
#include "./polynomial.h"


namespace uroboro {

	template<unsigned int N>
	void print_vec(vec<N> vec) {
		std::cout << "(";
		for (int i = 0; i < vec.size; ++i) {
			std::cout << vec[i];
			if(i != vec.size - 1)
				std::cout << ", ";
		}
		std::cout << ")";

		std::cout << std::endl;
	}

	template<unsigned int N, unsigned int K>
	void print_mat(mat<N, K> mat) {
		for (int i = 0; i < mat.row_size; ++i) {
				print_vec(mat.get_row(i));
		}

		std::cout << std::endl;
	}

	void print_complex(complex z) {

		std::cout << "(" << z.a;
		if(z.b == 1)
			std::cout << " + i)";
		else if(z.b >= 0)
			std::cout << " + " << z.b << "i)";
		else
			std::cout << " - " << abs(z.b) << "i)";

		std::cout << std::endl;
	}

	void print_complex_alg(complex z) {
		std::cout << "(" << z.a << ", " << z.b << ")" << std::endl;
	}

	void print_phasor(phasor z) {
		std::cout << z.modulus << "/" << z.phase << std::endl;
	}

	void print_quat(quat q) {
		std::cout << "(" << q.a << ", " << q.v[0] << ", " <<
				q.v[1] << ", " << q.v[2] << ")" << std::endl;
	}

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


	void print_polynomial(polynomial<real> p) {

		for (int i = 0; i < p.size(); ++i) {
			if(i) {
				std::cout << (p[i] >= 0 ? " + " : " - ") << uroboro::abs(p[i]) << "x^" << i;
			} else {
				std::cout << p[i];
			}	
		}

		std::cout << std::endl;
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
		std::cout << "Standard Deviation: " << sample_standard_deviation(X) << std::endl;
		std::cout << "Relative Error: " << sample_standard_relative_error(X) * 100 << "%" << std::endl;
		std::cout << "Mean Standard Deviation: " << sample_mean_standard_deviation(X) << std::endl;
		std::cout << "Chi-Square (Sigma): " << chi_square_sigma(X) << std::endl;
	}

	void print_population_stats(const vec_buff& X) {
		std::cout << "N = " << X.size() << std::endl;
		std::cout << "Mean: " << mean(X) << std::endl;
		std::cout << "Variance: " << variance(X) << std::endl;
		std::cout << "Standard Deviation: " << standard_deviation(X) << std::endl;
		std::cout << "Relative Error: " << standard_relative_error(X) * 100 << "%" << std::endl;
		std::cout << "Mean Standard Deviation: " << mean_standard_deviation(X) << std::endl;
		// std::cout << "Chi-Square (Sigma): " << chi_square_sigma(X) << std::endl;
	}

}

#endif
