#include "uroboro.h"

using namespace umath;


int main(int argc, char const *argv[]) {


	// Integral approximation
	std::cout << "Integral approximation" << std::endl;
	std::cout << approx_integral_midpoint(umath::square, 0, 3) << std::endl;
	std::cout << approx_integral_trapezoid(umath::square, 0, 3) << std::endl;
	std::cout << approx_integral_simpson(umath::square, 0, 3) << std::endl;

	std::cout << "Derivative approximation" << std::endl;
	std::cout << approx_derivative(square, 4) << std::endl;

	// Polynomials
	std::cout << "Polynomials\n" << std::endl;
	polynomial P1 = {1, 1, -1};
	polynomial P2 = {1, 2, 3};
	print_polynomial(P1);
	print_polynomial(P2);
	std::cout << std::endl;

	std::cout << "Approximated roots of P1 (different algorithms)" << std::endl;
	std::cout << approx_polyn_root_newton(P1, 1) << std::endl;
	std::cout << approx_polyn_root_steffensen(P1, 1) << std::endl;
	std::cout << approx_polyn_root_chebyshev(P1, 1) << std::endl;
	std::cout << std::endl;

	P1.trim();

	std::cout << "P1 + P2: ";
	print_polynomial(P1 + P2);
	std::cout << "P1 - P2: ";
	print_polynomial(P1 - P2);
	std::cout << "P1 * P2: ";
	print_polynomial(P1 * P2);
	std::cout << std::endl;

	std::cout << "d(P1)/dx: ";
	print_polynomial(differentiate_polynomial(P1));
	std::cout << "(P1 dx): ";
	print_polynomial(integrate_polynomial(P1));
	std::cout << std::endl;

	// Phasor
	std::cout << "Phasor\n" << std::endl;
	phasor p = phasor(PI, SQRT2);

	print_phasor(phasor(1, PI) * phasor(1, PI2));
	print_phasor(phasor(1, PI) / phasor(1, PI2));
	print_phasor(phasor(1, 0) + phasor(2, PI));
	std::cout << std::endl;

	// Quaternion rotation of a vector
	std::cout << "Quaternion rotation of a vector\n" << std::endl;
	vec3 v = {1, 2, 3};
	std::cout << "v = ";
	print_vec(v);

	std::cout << "quat::rotate(v, PI, {0, 1, 0}) = ";
	print_vec(quat::rotate(v, PI, {0, 1, 0}));
	std::cout << std::endl;

	// Complex product and rotation
	std::cout << "Complex product and rotation\n" << std::endl;
	complex z = {1, 0};
	complex w = {0, 1};
	complex r = complex::rotor(PI);

	std::cout << "z = ";
	print_complex(z);

	std::cout << "w = ";
	print_complex(w);

	std::cout << "r = ";
	print_complex(r);

	std::cout << "arg(r) = " << r.arg() << std::endl;

	std::cout << "z * w = ";
	print_complex(z * w);

	std::cout << "z * r = ";
	print_complex(z * r);
	std::cout << std::endl;

	// Complex functions
	std::cout << "Complex functions\n" << std::endl;
	std::cout << "sqrt(1, 1) = ";
	print_complex(sqrt({1, 1}));

	std::cout << "exp(1, 1) = ";
	print_complex(exp({1, 1}));

	std::cout << "sin(1, 1) = ";
	print_complex(sin({1, 1}));

	std::cout << "cos(1, 1) = ";
	print_complex(cos({1, 1}));

	std::cout << "tan(1, 1) = ";
	print_complex(tan({1, 1}));
	std::cout << std::endl;


	// Operations between vectors

	std::cout << "Operations between vectors\n" << std::endl;

	vec3 v1 = {1, 3, 1};
	vec3 v2 = {-1, 2, 1};

	std::cout << "v1 = ";
	print_vec(v1);

	std::cout << "v2 = ";
	print_vec(v2);

	std::cout << "v1 + v2 = ";
	print_vec(v1 + v2);

	std::cout << "v1 - v2 = ";
	print_vec(v1 - v2);

	std::cout << "v1 * v2 = " << (v1 * v2) << std::endl;

	std::cout << "v1 x v2 = ";
	print_vec(v1.cross(v2));
	std::cout << std::endl;

	// Matrix operations
	std::cout << "Matrix operations\n" << std::endl;

	mat4 m1 = mat4(2);
	m1.at(1, 2) = 3;
	m1.at(2, 1) = 5;

	mat<3, 4> m2 = mat<3, 4>(3);

	print_mat(m1);
	print_mat(m1.transposed());
	print_mat(m1 * m2);

	print_mat(mat4::identity() * mat4::diagonal(2.0));
	print_mat(mat4::translation(vec3({1, 2, 3})));


	// Matrix rotation around an arbitrary axis
	v = {1, 1, 1};
	mat3 R = mat3::rotation_3x3(PI2, vec3({1, 1, 1}));

	print_vec(R * v);

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
