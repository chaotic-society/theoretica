#ifndef UROBORO_UTILITY_H
#define UROBORO_UTILITY_H
#include <iostream>

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

}

#endif
