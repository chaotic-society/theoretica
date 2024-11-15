
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


constexpr unsigned int DEFAULT_ITER = 10;


// Compute the L_inf norm of any iterable structure, such as vectors or matrices.
// The norm finds the maximum element in absolute value.
template<typename Structure>
real linf_norm(const Structure& A) {

	real m = 0.0;

	for (const auto& x : A)
		m = std::max(m, std::abs(x));

	return m;
}


// Generate a random vector with Gaussian distributed elements
template<typename Vector = vec<real>>
Vector rand_vec(real m, real s, unsigned int n) {

	Vector v;
	v.resize(n);

	for (auto& x : v)
		x = random::gaussian(m, s);

	return v;
}


// Generate a random matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix A;
	A.resize(rows, cols);

	for (auto& x : A)
		x = random::gaussian(m, s);

	return A;
}


// Generate a random lower triangular matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat_lower(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix L;
	L.resize(rows, cols);

	for (unsigned int i = 0; i < L.rows(); ++i)
		for (unsigned int j = 0; j < L.cols(); ++j)
			if (i <= j)
				L(i, j) = random::gaussian(m, s);

	return L;
}


// Generate a random upper triangular matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat_upper(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix U;
	U.resize(rows, cols);

	for (unsigned int i = 0; i < U.rows(); ++i)
		for (unsigned int j = 0; j < U.cols(); ++j)
			if (i >= j)
				U(i, j) = random::gaussian(m, s);

	return U;
}


// Generate a random symmetric matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat_symmetric(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix A = rand_mat<Matrix>(m, s, rows, cols);
	return (A + algebra::transpose(A)) * 0.5;
}


// Generate a random positive definite symmetric matrix with random elements
template<typename Matrix = mat<real>>
Matrix rand_mat_posdef(real m, real s, unsigned int rows) {

	Matrix A = rand_mat<Matrix>(m, s, rows, rows);
	return algebra::mat_mul_transpose(A, A);
}


// Estimate error of a function over matrices. The residual function is expected to
// generate a casual matrix and compute the residual of the tested function, while
// the expected function returns the expected result of the computation (usually zero).
template<typename Matrix = mat<real>>
auto mat_estimator() {

	return [](
		std::function<real()> residual,
		std::function<real()> expected,
		prec::estimate_options<real> opt) -> prec::estimate_result {

		long double absErr = 0.0;
		long double sqrAbsErr = 0.0;
		long double maxErr = -inf();

		for (unsigned int i = 0; i < opt.iterations; ++i) {

			const long double r = std::abs(residual() - expected());
			
			absErr += r;
			sqrAbsErr += r * r;
			maxErr = std::max(maxErr, r);
		}

		prec::estimate_result res;
		res.maxErr = maxErr;
		res.meanErr = absErr / opt.iterations;
		res.rmsErr = std::sqrt(sqrAbsErr / opt.iterations);
		return res;
	};
}


// Run a test against a residual function over random matrices or vectors.
template<typename Function>
void test_residual(
	const std::string& name,
	Function residual,
	unsigned int iterations = DEFAULT_ITER) {

	auto opt = prec::estimate_options<real>(
		prec::interval(),
		mat_estimator()
	);
	opt.iterations = iterations;

	prec::estimate(
		name,
		residual,
		[]() { return 0.0; },
		opt
	);
}


int main(int argc, char const *argv[]) {

	prec::settings.outputFiles = { "test/prec/test_algebra.csv" };

	prec::setup("algebra", argc, argv);

		// algebra.h


		const unsigned int N = 100;


	{
		vec<real> v = vec<real>(N);
		algebra::vec_error(v);

		prec::equals(
			"vec_error",
			is_nan(v),
			true, 0
		);
	}


	{
		mat<real> A = mat<real>(N, N);
		algebra::mat_error(A);

		prec::equals(
			"mat_error",
			is_nan(A),
			true, 0
		);
	}


		test_residual("normalize", []() {

			auto v = rand_vec(0.0, 1.0, N);
			auto w = algebra::normalize(v);
			return std::abs(1 - algebra::norm(w));
		});


		test_residual("make_normalized", []() {

			auto v = rand_vec(0.0, 1.0, N);
			algebra::make_normalized(v);
			return std::abs(1 - algebra::norm(v));
		});


		test_residual("dot", []() {

			auto v = rand_vec(0.0, 1.0, N);
			return std::abs(algebra::dot(v, v) - algebra::sqr_norm(v));
		});


		test_residual("cross", []() {
			auto v1 = rand_vec(0.0, 1.0, 3);
			auto v2 = rand_vec(0.0, 1.0, 3);
			return std::abs(v1 * algebra::cross(v1, v2));
		});


		test_residual("cross", []() {
			auto v1 = rand_vec(0.0, 1.0, 3);
			auto v2 = rand_vec(0.0, 1.0, 3);
			return std::abs(v2 * algebra::cross(v1, v2));
		});


		test_residual("transpose", []() {

			auto A = rand_mat(0.0, 1.0, N, N);
			return linf_norm(A - algebra::transpose(algebra::transpose(A)));
		});


		test_residual("make_transposed", []() {

			auto A = rand_mat(0.0, 1.0, N, N);
			auto B = A;

			algebra::make_transposed(B);
			algebra::make_transposed(B);
			return linf_norm(A - B);
		});


		test_residual("decompose_cholesky", []() {

			auto A = rand_mat_posdef(0.0, 1.0, N);
			auto L = algebra::decompose_cholesky(A);
			return linf_norm(A - algebra::mat_mul_transpose(L, L));
		});


		test_residual("decompose_cholesky_inplace", []() {

			auto A = rand_mat_posdef(0.0, 1.0, N);
			auto L = A;
		
			algebra::decompose_cholesky_inplace(L);
			return linf_norm(A - algebra::mat_mul_transpose(L, L));
		});


		test_residual("det", []() {

			size_t sz = 10;

			auto L = rand_mat_lower(0.0, 1.0, sz, sz);
			auto U = rand_mat_upper(0.0, 1.0, sz, sz);
			auto A = L * U;

			matrix_element_t<decltype(A)> d = 1.0;

			// The determinant of L * U is the product
			// of the diagonal elements of L and U
			for (size_t i = 0; i < sz; ++i)
				d *= L(i, i) * U(i, i);

			return std::abs(d - algebra::det(A));
		});


		// mat.h

		// vec.h

		// distance.h
	
	prec::terminate();
}
