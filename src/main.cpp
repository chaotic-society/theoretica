#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace umath;

int main(int argc, char const *argv[]) {

	complex z = complex(1, 0);
	complex w = complex(0, 1);
	print_complex_alg(z);
	print_complex_alg(w);
	print_complex_alg(z * w);
	print_complex_alg(z * w * w);
	print_complex_alg(z / w);
	print_complex_alg(z / w / w);

	auto vec1 = vec4({1, 3, 1, 2});
	auto vec2 = vec4({-1, 2, 1, 0});

	auto mat1 = mat4(2);
	mat1.at(1, 2) = 3;
	mat1.at(2, 1) = 5;
	print_mat(mat1);
	print_mat(mat1.transposed());

	// std::cout << mat1.is_square() << std::endl;
	// std::cout << mat1.is_symmetric() << std::endl;
	// std::cout << mat1.is_diagonal() << std::endl;
	// mat1.transpose();
	// print_mat(mat1);


	// auto X = vec<3>({1, 1.2, 0.9});
	// auto Y = vec<3>({10.1, 11.7, 9.2});

	// std::cout << covariance(X, Y) << std::endl;
	// std::cout << sample_standard_deviation(X) << std::endl;
	// std::cout << sample_standard_deviation(Y) << std::endl;
	// std::cout << mean_standard_deviation(X) << std::endl;
	// std::cout << mean_standard_deviation(Y) << std::endl;

	return 0;
}
