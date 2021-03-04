#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace umath;

int main(int argc, char const *argv[]) {

	auto vec1 = vec<4>({1, 3, 1, 2});
	auto vec2 = vec<4>({-1, 2, 1, 0});

	auto mat1 = mat<4, 3>(2);
	auto mat2 = mat1 * 3;
	print_vec(mat1 * vec1);
	print_vec(mat2 * vec2);

	std::cout << mat1.is_square() << std::endl;
	std::cout << mat1.is_symmetric() << std::endl;
	std::cout << mat1.is_diagonal() << std::endl;

	auto X = vec<3>({1, 1.2, 0.9});
	auto Y = vec<3>({10.1, 11.7, 9.2});

	std::cout << covariance(X, Y) << std::endl;
	std::cout << sample_standard_deviation(X) << std::endl;
	std::cout << sample_standard_deviation(Y) << std::endl;
	std::cout << mean_standard_deviation(X) << std::endl;
	std::cout << mean_standard_deviation(Y) << std::endl;

	return 0;
}
