#ifndef UROBORO_UTILITY_H
#define UROBORO_UTILITY_H
#include <iostream>

#include "./vec.h"
#include "./mat.h"
#include "./complex.h"

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
			std::cout << " + " << z.b << ")";
		else
			std::cout << " - " << z.b << ")";

		std::cout << std::endl;
	}

	void print_complex_alg(complex z) {
		std::cout << "(" << z.a << ", " << z.b << ")" << std::endl;
	}

}

#endif
