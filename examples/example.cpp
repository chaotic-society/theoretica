
///
/// @file example.cpp Example file
///

#include "uroboro.h"

using namespace umath;


int main(int argc, char const *argv[]) {

	std::cout.precision(12);

	// Automatic differentiation using dual numbers
	std::cout << "Automatic differentiation" << std::endl;

	dual d = dual(2, 1);
	std::cout << square(d) << std::endl;

	// Integral approximation
	std::cout << "Integral approximation" << std::endl;
	std::cout << approx_integral_midpoint(square, 0, 3) << std::endl;
	std::cout << approx_integral_trapezoid(square, 0, 3) << std::endl;
	std::cout << approx_integral_simpson(square, 0, 3) << std::endl;
	std::cout << std::endl;

	std::cout << "Derivative approximation" << std::endl;
	std::cout << approx_derivative_central(square, 4) << std::endl;
	std::cout << approx_derivative_forward(square, 4) << std::endl;
	std::cout << std::endl;

	// Polynomials
	std::cout << "Polynomials\n" << std::endl;
	polynomial<real> P1 = {1, 1, -1};
	polynomial<real> P2 = {1, 2, 3};
	std::cout << P1 << std::endl << P2 << std::endl;

	std::cout << "Approximated roots of P1 (different algorithms)" << std::endl;
	std::cout << approx_polyn_root_newton(P1, 1) << std::endl;
	std::cout << approx_polyn_root_steffensen(P1, 1) << std::endl;
	std::cout << approx_polyn_root_chebyshev(P1, 1) << std::endl;
	std::cout << std::endl;

	P1.trim();

	std::cout << "P1 + P2: ";
	std::cout << (P1 + P2) << std::endl;
	std::cout << "P1 - P2: ";
	std::cout << (P1 - P2) << std::endl;
	std::cout << "P1 * P2: ";
	std::cout << (P1 * P2) << std::endl;
	std::cout << std::endl;

	std::cout << "d(P1)/dx: ";
	std::cout << differentiate_polynomial(P1) << std::endl;
	std::cout << "(P1 dx): ";
	std::cout << integrate_polynomial(P1) << std::endl;
	std::cout << std::endl;

	// Phasor
	std::cout << "Phasor\n" << std::endl;
	phasor p = phasor(PI, SQRT2);

	std::cout << phasor(1, PI) * phasor(1, PI2) << std::endl;
	std::cout << phasor(1, PI) / phasor(1, PI2) << std::endl;
	std::cout << phasor(1, 0) + phasor(2, PI) << std::endl;
	std::cout << std::endl;

	// Quaternion rotation of a vector
	std::cout << "Quaternion rotation of a vector\n" << std::endl;
	vec3 v = {1, 2, 3};
	std::cout << "v = " << v;

	std::cout << "quat::rotate(v, PI, {0, 1, 0}) = "
		<< quat::rotate(v, PI, {0, 1, 0}) << std::endl;

	// Complex product and rotation
	std::cout << "Complex product and rotation\n" << std::endl;
	complex z = {1, 0};
	complex w = {0, 1};
	complex r = complex::rotor(PI);

	std::cout << "z = " << z << std::endl;
	std::cout << "w = " << w << std::endl;
	std::cout << "r = " << r << std::endl;
	std::cout << "arg(r) = " << r.arg() << std::endl;
	std::cout << "z * w = " << (z * w) << std::endl;
	std::cout << "z * r = " << (z * r) << std::endl;
	std::cout << std::endl;

	// Complex functions
	std::cout << "Complex functions\n" << std::endl;
	std::cout << "sqrt(1, 1) = ";
	std::cout << sqrt(z + w) << std::endl;

	std::cout << "exp(1, 1) = ";
	std::cout << exp(z + w) << std::endl;

	std::cout << "sin(1, 1) = ";
	std::cout << sin(z + w) << std::endl;

	std::cout << "cos(1, 1) = ";
	std::cout << cos(z + w) << std::endl;

	std::cout << "tan(1, 1) = ";
	std::cout << tan(z + w) << std::endl;
	std::cout << std::endl;


	// Operations between vectors

	std::cout << "Operations between vectors\n" << std::endl;

	vec3 v1 = {1, 3, 1};
	vec3 v2 = {-1, 2, 1};

	std::cout << "v1 = " << v1 << std::endl;
	std::cout << "v2 = " << v2 << std::endl;
	std::cout << "v1 + v2 = " << (v1 + v2) << std::endl;
	std::cout << "v1 - v2 = " << (v1 - v2) << std::endl;
	std::cout << "v1 * v2 = " << (v1 * v2) << std::endl;
	std::cout << "v1 x v2 = " << v1.cross(v2) << std::endl;
	std::cout << std::endl;

	// Matrix operations
	std::cout << "Matrix operations\n" << std::endl;

	mat4 m1 = mat4(2);
	m1.at(1, 2) = 3;
	m1.at(2, 1) = 5;

	mat<4, 3> m2 = mat<4, 3>(3);

	std::cout << m1 << std::endl;
	std::cout << (m1.transposed()) << std::endl;
	std::cout << (m1 * m2) << std::endl;

	std::cout << (mat4::identity() * mat4::diagonal(2.0)) << std::endl;
	std::cout << mat4::translation(vec3({1, 2, 3})) << std::endl;


	// Matrix rotation around an arbitrary axis
	v = {1, 1, 1};
	mat3 R = mat3::rotation_3x3(PI2, vec3({1, 1, 1}));

	std::cout << (R * v) << std::endl;
	std::cout << std::endl;


	// Statistical analysis
	std::cout << "Statistical analysis\n" << std::endl;

	vec_buff X = {1, 1.2, 0.9, 0.78, 0.71};
	vec_buff Y = {10.1, 11.7, 9.2, 8.1, 6.9};

	std::cout << "X = ";
	print_vec_buff_row(X);

	std::cout << "Y = ";
	print_vec_buff_row(Y);

	std::cout << std::endl;

	std::cout << "covar(X, Y) = " << sample_covariance(X, Y) << std::endl;
	std::cout << "correlation(X, Y) = " << sample_correlation_coefficient(X, Y) << std::endl;
	std::cout << std::endl;

	std::cout << "stdev(X) = " << smpl_stdev(X) << std::endl;
	std::cout << "relative_err(X) = " << sample_standard_relative_error(X) << std::endl;
	std::cout << std::endl;

	std::cout << "stdev(Y) = " << smpl_stdev(Y) << std::endl;
	std::cout << "relative_err(Y) = " << sample_standard_relative_error(Y) << std::endl;
	std::cout << std::endl;

	std::cout << "STDOM(X) = " << smpl_stdom(X) << std::endl;
	std::cout << "STDOM(Y) = " << smpl_stdom(Y) << std::endl;
	std::cout << std::endl;

	std::cout << lst_sqrs_lin_intercept(X, Y) << std::endl;
	std::cout << lst_sqrs_lin_slope(X, Y) << std::endl;
	std::cout << std::endl;

	std::cout << distribution::gaussian(1.0, 1.0, 0.02) << std::endl;
	std::cout << distribution::binomial(0, 4, 0.5) << std::endl;

	return 0;
}
