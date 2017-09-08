#ifndef UROBORO_UTILITY_H
#define UROBORO_UTILITY_H
#include <iostream>
#include <ctime>
#include "mat.h"
#include "vec.h"


namespace uroboro {

		inline void set_output_prec(unsigned int prec) {
			std::cout.precision(prec);
		}

		inline void print_mat(mat4 matrix) {

			std::cout << matrix.data[0][0] << ", " << matrix.data[1][0] << ", " <<
			matrix.data[2][0] << ", " << matrix.data[3][0] << std::endl;

			std::cout << matrix.data[0][1] << ", " << matrix.data[1][1] << ", " <<
			matrix.data[2][1] << ", " << matrix.data[3][1] << std::endl;

			std::cout << matrix.data[0][2] << ", " << matrix.data[1][2] << ", " <<
			matrix.data[2][2] << ", " << matrix.data[3][2] << std::endl;

			std::cout << matrix.data[0][3] << ", " << matrix.data[1][3] << ", " <<
			matrix.data[2][3] << ", " << matrix.data[3][3] << std::endl;
		}

		inline void print_vec(vec4 vector) {
			std::cout << vector.x << ", " << vector.y << ", " << vector.z << ", " << vector.w << std::endl;
		}

		inline void start_clock(std::clock_t& clock) {
			clock = std::clock();
		}

		inline unsigned int end_clock(std::clock_t& begin) {
			std::clock_t end = std::clock();
			return double(end - begin) / CLOCKS_PER_SEC;
		}

}

#endif
