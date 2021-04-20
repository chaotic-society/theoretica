#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace umath;

int main(int argc, char const *argv[]) {

	// Quaternion rotation of a vector
	vec3 v = {1, 2, 3};
	print_vec(v);
	print_vec(quat::rotate(v, PI, {0, 1, 0}));

	// Complex multiplication and rotation
	complex z = {1, 0};
	complex w = {0, 1};
	complex r = complex::rotor(PI);

	print_complex(z);
	print_complex(w);
	print_complex(z * w);
	print_complex(z * r);

	// Dot product and cross product between vectors
	vec3 v1 = {1, 3, 1};
	vec3 v2 = {-1, 2, 1};

	std::cout << v1 * v2 << std::endl;
	print_vec(v1.cross(v2));

	// Matrix product and transposition
	auto m1 = mat4(2);
	auto m2 = mat<3, 4>(3);
	m1.at(1, 2) = 3;
	m1.at(2, 1) = 5;

	print_mat(m1);
	print_mat(m1.transposed());
	print_mat(m1 * m2);

	print_mat(mat4::identity() * mat4::diagonal(2.0));

	// Statistical analysis
	vec_buff X = {1, 1.2, 0.9, 0.78, 0.71};
	vec_buff Y = {10.1, 11.7, 9.2, 8.1, 6.9};

	std::cout << sample_covariance(X, Y) << std::endl;
	std::cout << smpl_stdev(X) << std::endl;
	std::cout << sample_standard_relative_error(X) << std::endl;
	std::cout << smpl_stdev(Y) << std::endl;
	std::cout << sample_standard_relative_error(Y) << std::endl;
	std::cout << smpl_stdom(X) << std::endl;
	std::cout << smpl_stdom(Y) << std::endl;
	std::cout << lst_sqrs_lin_intercept(X, Y) << std::endl;
	std::cout << lst_sqrs_lin_slope(X, Y) << std::endl;

	return 0;
}
