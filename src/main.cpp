#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace umath;

int main(int argc, char const *argv[]) {

	vec3 v = {1, 2, 3};
	print_vec(v);
	print_vec(quat::rotate(v, PI, {0, 1, 0}));

	quat q = {1, 0, 0, 1};

	complex z = {1, 0};
	complex w = {0, 1};
	print_complex_alg(z);
	print_complex_alg(w);
	print_complex_alg(z * w);

	vec4 vec1 = {1, 3, 1, 2};
	vec4 vec2 = {-1, 2, 1, 0};

	auto m1 = mat4(2);
	auto m2 = mat<3, 4>(3);
	m1.at(1, 2) = 3;
	m1.at(2, 1) = 5;
	print_mat(m1);
	print_mat(m1.transposed());
	print_mat(m1 * m2);

	print_mat(mat4::identity() * mat4::diagonal(2.0));

	vec_buff X = {1, 1.2, 0.9};
	vec_buff Y = {10.1, 11.7, 9.2};

	std::cout << covariance(X, Y) << std::endl;
	std::cout << sample_standard_deviation(X) << std::endl;
	std::cout << sample_standard_deviation(Y) << std::endl;
	std::cout << mean_standard_deviation(X) << std::endl;
	std::cout << mean_standard_deviation(Y) << std::endl;

	return 0;
}
