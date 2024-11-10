
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


// Compute the L_inf norm of any iterable structure, such as vectors or matrices.
// The norm finds the maximum element in absolute value.
template<typename Structure>
real linf_norm(const Structure& A) {

	real m = 0.0;

	for (const auto& x : A)
		m = std::max(m, std::abs(x));

	return m;
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


int main(int argc, char const *argv[]) {

	prec::settings.outputFiles = { "test/prec/test_algebra.csv" };

	prec::setup("algebra", argc, argv);

		// algebra.h

		auto alg_opt = prec::estimate_options<real>(
			prec::interval(),
			mat_estimator()
		);

		// Compute the residual N items
		alg_opt.iterations = 100;

		auto is_zero = []() -> real { return 0.0; };


		prec::estimate(
			"transpose",
			[]() {
				auto A = rand_mat(0.0, 1.0, 100, 100);
				return linf_norm(A - algebra::transpose(algebra::transpose(A)));
			},
			is_zero, alg_opt
		);


		prec::estimate(
			"make_transposed",
			[]() {
				auto A = rand_mat(0.0, 1.0, 100, 100);
				auto B = A;
				algebra::make_transposed(B);
				algebra::make_transposed(B);
				return linf_norm(A - B);
			},
			is_zero, alg_opt
		);


		prec::estimate(
			"decompose_cholesky",
			[]() {
				auto A = rand_mat_posdef(0.0, 1.0, 100);
				auto L = algebra::decompose_cholesky(A);
				return linf_norm(A - algebra::mat_mul_transpose(L, L));
			},
			is_zero, alg_opt
		);


		prec::estimate(
			"decompose_cholesky_inplace",
			[]() {
				auto A = rand_mat_posdef(0.0, 1.0, 100);
				auto L = A;
				algebra::decompose_cholesky_inplace(L);
				return linf_norm(A - algebra::mat_mul_transpose(L, L));
			},
			is_zero, alg_opt
		);

		// mat.h

		// vec.h

		// distance.h
	
	prec::terminate();
}
