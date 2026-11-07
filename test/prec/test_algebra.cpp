
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


constexpr unsigned int DEFAULT_ITER = 1000;


// Source of random numbers
random::random_source rnd {0};


// Compute the L_inf norm of any iterable structure, such as vectors or matrices.
// The norm finds the maximum element in absolute value.
template<typename Structure>
real linf_norm(const Structure& A) {

	real m = 0.0;

	for (const auto& x : A) {

		const auto abs_x = std::abs(x);

		if (abs_x != abs_x)
			return inf();

		m = std::max(m, abs_x);
	}

	return m;
}


// Generate a random vector with Gaussian distributed elements
template<typename Vector = vec<real>>
Vector rand_vec(real m, real s, unsigned int n) {

	Vector v;
	v.resize(n);

	for (auto& x : v)
		x = rnd.gaussian(m, s);

	return v;
}


// Generate a random matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix A;
	A.resize(rows, cols);

	for (auto& x : A)
		x = rnd.gaussian(m, s);

	return A;
}


// Generate a random lower triangular matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat_lower(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix L;
	L.resize(rows, cols);

	for (unsigned int i = 0; i < L.rows(); ++i)
		for (unsigned int j = 0; j < L.cols(); ++j)
			if (i >= j)
				L(i, j) = rnd.gaussian(m, s);

	return L;
}


// Generate a random upper triangular matrix with Gaussian distributed elements
template<typename Matrix = mat<real>>
Matrix rand_mat_upper(real m, real s, unsigned int rows, unsigned int cols) {

	Matrix U;
	U.resize(rows, cols);

	for (unsigned int i = 0; i < U.rows(); ++i)
		for (unsigned int j = 0; j < U.cols(); ++j)
			if (i <= j)
				U(i, j) = rnd.gaussian(m, s);

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


// Custom matrix-matrix multiplication
template<typename Matrix1, typename Matrix2>
auto mat_mul(const Matrix1& A, const Matrix2& B) {

	using Result = decltype(A * B);
	Result C;
	C.resize(A.rows(), B.cols());
	algebra::mat_zeroes(C);

	for (unsigned int i = 0; i < A.rows(); ++i)
		for (unsigned int j = 0; j < B.cols(); ++j)
			for (unsigned int k = 0; k < A.cols(); ++k)
				C(i, j) += A(i, k) * B(k, j);

	return C;
}


// Custom vector-matrix multiplication
template<typename Vector, typename Matrix>
auto vec_mat_mul(const Vector& v, const Matrix& A) {

	using Result = decltype(v * A);
	Result r;
	r.resize(A.cols());
	algebra::vec_zeroes(r);

	for (unsigned int j = 0; j < A.cols(); ++j)
		for (unsigned int i = 0; i < A.rows(); ++i)
			r[j] += v[i] * A(i, j);

	return r;
}


// Custom matrix-vector multiplication
template<typename Matrix, typename Vector>
auto mat_vec_mul(const Matrix& A, const Vector& v) {

	using Result = decltype(A * v);
	Result r;
	r.resize(A.rows());
	algebra::vec_zeroes(r);

	for (unsigned int i = 0; i < A.rows(); ++i)
		for (unsigned int j = 0; j < A.cols(); ++j)
			r[i] += A(i, j) * v[j];

	return r;
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
	prec::prec_context& ctx,
	const std::string& name,
	Function residual,
	unsigned int iterations = DEFAULT_ITER) {

	auto opt = prec::estimate_options<real>(
		prec::interval(),
		mat_estimator()
	);
	opt.iterations = iterations;

	ctx.estimate(
		name,
		residual,
		[]() { return 0.0; },
		opt
	);
}


int main(int argc, char const *argv[]) {

	auto ctx = prec::make_context("algebra", argc, argv);
	ctx.settings.outputFiles = { "test/prec/test_algebra.csv" };
	ctx.settings.multithreading = false;
	rnd = ctx.random->get_rnd();


	// Equation options for vector comparison
	prec::equation_options<vec4> vec4_opt { 1E-08, prec::distance::euclidean<vec4> };
	prec::equation_options<vec3> vec3_opt { 1E-08, prec::distance::euclidean<vec3> };
	
	
	// algebra.h

	const unsigned int N = 16;


	{
		vec<real> v = vec<real>(N);
		algebra::vec_error(v);

		ctx.equals(
			"vec_error",
			is_nan(v),
			true, 0
		);
	}


	{
		mat<real> A = mat<real>(N, N);
		algebra::mat_error(A);

		ctx.equals(
			"mat_error",
			is_nan(A),
			true, 0
		);
	}


	test_residual(ctx, "normalize", []() {

		auto v = rand_vec(0.0, 1.0, N);
		auto w = algebra::normalize(v);
		return std::abs(1 - algebra::norm(w));
	});


	test_residual(ctx, "make_normalized", []() {

		auto v = rand_vec(0.0, 1.0, N);
		algebra::make_normalized(v);
		return std::abs(1 - algebra::norm(v));
	});


	test_residual(ctx, "dot", []() {

		auto v = rand_vec(0.0, 1.0, N);
		return std::abs(algebra::dot(v, v) - algebra::sqr_norm(v));
	});


	test_residual(ctx, "cross", []() {
		auto v1 = rand_vec(0.0, 1.0, 3);
		auto v2 = rand_vec(0.0, 1.0, 3);
		return std::abs(v1 * algebra::cross(v1, v2));
	});


	test_residual(ctx, "cross", []() {
		auto v1 = rand_vec(0.0, 1.0, 3);
		auto v2 = rand_vec(0.0, 1.0, 3);
		return std::abs(v2 * algebra::cross(v1, v2));
	});


	test_residual(ctx, "transpose", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		return linf_norm(A - algebra::transpose(algebra::transpose(A)));
	});


	test_residual(ctx, "make_transposed", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto B = A;

		algebra::make_transposed(B);
		algebra::make_transposed(B);
		return linf_norm(A - B);
	});


	test_residual(ctx, "decompose_cholesky", []() {

		auto A = rand_mat_posdef(0.0, 1.0, N);
		auto L = algebra::decompose_cholesky(A);
		return linf_norm(A - algebra::mat_mul_transpose(L, L));
	});


	test_residual(ctx, "decompose_cholesky_inplace", []() {

		auto A = rand_mat_posdef(0.0, 1.0, N);
		auto L = A;
	
		algebra::decompose_cholesky_inplace(L);
		return linf_norm(A - algebra::mat_mul_transpose(L, L));
	});


	test_residual(ctx, "det", []() {

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


	test_residual(ctx, "mat_sum", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto B = rand_mat(0.0, 1.0, N, N);
		
		mat<real> C;
		C.resize(N, N);
		algebra::mat_sum(C, A, B);

		return linf_norm(C - A - B);
	});


	test_residual(ctx, "mat_diff", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto B = rand_mat(0.0, 1.0, N, N);
		mat<real> C;
		C.resize(N, N);
		algebra::mat_diff(C, A, B);

		return linf_norm(C - A + B);
	});


	test_residual(ctx, "mat_lincomb", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto B = rand_mat(0.0, 1.0, N, N);

		const real a = 1.5;
		const real b = -0.75;

		mat<real> C;
		C.resize(N, N);
		algebra::mat_lincomb(C, a, A, b, B);

		mat<real> expected;
		expected.resize(N, N);
		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < N; ++j)
				expected(i, j) = a * A(i, j) + b * B(i, j);

		return linf_norm(C - expected);
	});


	test_residual(ctx, "mat_mul", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto B = rand_mat(0.0, 1.0, N, N);
		auto C = algebra::mat_mul(A, B);

		return linf_norm(C - mat_mul(A, B));
	});


	test_residual(ctx, "vec_mat_mul", []() {

		auto v = rand_vec(0.0, 1.0, N);
		auto A = rand_mat(0.0, 1.0, N, N);
		auto r = algebra::vec_mat_mul<vec<real>>(v, A);

		return linf_norm(r - vec_mat_mul(v, A));
	});


	test_residual(ctx, "apply_transform", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto v = rand_vec(0.0, 1.0, N);
		auto r = v;
		algebra::apply_transform(A, r);

		return linf_norm(r - mat_vec_mul(A, v));
	});


	test_residual(ctx, "mat_mul_transpose", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto R = algebra::mat_mul_transpose(A, A);

		return linf_norm(R - algebra::transpose(R));
	});


	test_residual(ctx, "mat_transpose_mul", []() {

		auto A = rand_mat(0.0, 1.0, N, N);
		auto R = algebra::mat_transpose_mul(A, A);
		return linf_norm(R - algebra::transpose(R));
	});


	test_residual(ctx, "decompose_lu", []() {

		auto L = rand_mat_lower(0.0, 1.0, N, N);
		auto U = rand_mat_upper(0.0, 1.0, N, N);

		// Diagonal is assumed to be 1 for L
		for (unsigned int i = 0; i < N; ++i)
			L(i, i) = 1.0;

		auto A = L * U;
		mat<real> L2;
		mat<real> U2;
		algebra::decompose_lu(A, L2, U2);

		return linf_norm(A - mat_mul(L2, U2));
	});


	test_residual(ctx, "decompose_lu_inplace", []() {

		auto L = rand_mat_lower(0.0, 1.0, N, N);
		auto U = rand_mat_upper(0.0, 1.0, N, N);

		for (unsigned int i = 0; i < N; ++i)
			L(i, i) = 1.0;

		auto A = L * U;
		auto B = A;
		algebra::decompose_lu_inplace(B);

		mat<real> L2 (N, N);
		mat<real> U2 (N, N);

		algebra::mat_zeroes(L2);
		algebra::mat_zeroes(U2);

		for (unsigned int i = 0; i < N; ++i) {

			L2(i, i) = 1.0;

			for (unsigned int j = 0; j < i; ++j)
				L2(i, j) = B(i, j);

			for (unsigned int j = i; j < N; ++j)
				U2(i, j) = B(i, j);
		}

		return linf_norm(A - mat_mul(L2, U2));
	});


	test_residual(ctx, "solve_triangular_lower", []() {

		auto L = rand_mat_lower(0.0, 1.0, N, N);
		for (unsigned int i = 0; i < N; ++i)
			L(i, i) += 5.0;

		auto x_true = rand_vec(0.0, 1.0, N);
		auto b = mat_vec_mul(L, x_true);
		auto x = algebra::solve_triangular_lower(L, b);

		return linf_norm(x - x_true);
	});


	test_residual(ctx, "solve_triangular_upper", []() {

		auto U = rand_mat_upper(0.0, 1.0, N, N);
		for (unsigned int i = 0; i < N; ++i)
			U(i, i) += 5.0;

		auto x_true = rand_vec(0.0, 1.0, N);
		auto b = mat_vec_mul(U, x_true);
		auto x = algebra::solve_triangular_upper(U, b);

		return linf_norm(x - x_true);
	});


	test_residual(ctx, "solve_lu", []() {

		auto L = rand_mat_lower(0.0, 1.0, N, N);
		auto U = rand_mat_upper(0.0, 1.0, N, N);

		for (unsigned int i = 0; i < N; ++i) {
			L(i, i) = 1.0;
			U(i, i) += 5.0;
		}
		
		auto A = L * U;
		auto x_true = rand_vec(0.0, 1.0, N);
		auto b = mat_vec_mul(A, x_true);
		auto x = algebra::solve_lu(A, b);
		
		return linf_norm(x - x_true);
	});


	test_residual(ctx, "solve_cholesky", []() {

		auto A = rand_mat_posdef(0.0, 1.0, N);

		// Make diagonally dominant
		for (unsigned int i = 0; i < N; ++i)
			A(i, i) += 5.0;

		auto L = algebra::decompose_cholesky(A);
		auto x_true = rand_vec(0.0, 1.0, N);
		auto b = mat_vec_mul(A, x_true);
		auto x = algebra::solve_cholesky(L, b);

		return linf_norm(x - x_true);
	});


	test_residual(ctx, "inverse", []() {

		auto L = rand_mat_lower(0.0, 1.0, N, N);
		auto U = rand_mat_upper(0.0, 1.0, N, N);

		for (unsigned int i = 0; i < N; ++i) {
			L(i, i) = 1.0;
			U(i, i) += 5.0;
		}

		auto A = L * U;
		auto inv = algebra::inverse(A);
		auto I = algebra::identity<mat<real>>(N, N);

		return linf_norm(mat_mul(A, inv) - I);
	});


	{
		mat4 A = algebra::identity<mat4>();
		algebra::mat_shift_diagonal(A, 2.5);

		ctx.equals("mat_shift_diagonal", A(0, 0), 3.5);
		ctx.equals("mat_shift_diagonal", A(1, 1), 3.5);
		ctx.equals("mat_shift_diagonal", A(0, 1), 0.0);
	}


	// TO-DO: Implement random orthogonal matrices for testing eigensolvers

	{

		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {0.0, 0.0, 1.0};

		auto lambda = algebra::eigenvalue_power(A, x, 1E-08);

		ctx.equals("eigenvalue_power", lambda, 5.0, 1E-08);
	}


	{

		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {0.0, 0.0, 1.0};
		vec3 v;
		
		auto lambda = algebra::eigenpair_power(A, x, v, 1E-08);
		ctx.equals("eigenpair_power", mat_vec_mul(A, v), v * lambda, vec3_opt);
	}


	{

		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {1.0, 0.0, 0.0};
		auto lambda = algebra::eigenvalue_inverse(A, x, 1E-08);
		
		ctx.equals("eigenvalue_inverse", lambda, 1.0, 1E-08);
	}


	{

		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {1.0, 0.0, 0.0};
		vec3 v;

		auto lambda = algebra::eigenpair_inverse(A, x, v, 1E-08);
		
		ctx.equals("eigenpair_inverse", mat_vec_mul(A, v), v * lambda, vec3_opt);
	}


	{

		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {0.0, 1.0, 0.0};
		
		auto lambda = algebra::eigenvalue_rayleigh(A, 1.2, x, 1E-08, 1000);
		
		ctx.equals("eigenvalue_rayleigh", lambda, 2.0, 1E-08);
	}


	{
		mat3 A = {
			{1.0, 0.0, 0.0},
			{0.0, 2.0, 0.0},
			{0.0, 0.0, 5.0}
		};

		vec3 x = {0.0, 1.0, 0.0};

		ctx.equals("rayleigh_quotient", algebra::rayleigh_quotient(A, x), 2.0, 1E-08);
	}


	// algebra_types.h
	{

		ctx.equals("vec2::size()", vec2(2).size(), 2);
		ctx.equals("vec3::size()", vec3(3).size(), 3);
		ctx.equals("vec4::size()", vec4(4).size(), 4);

		ctx.equals("mat2::rows()", mat2(0.0).rows(), 2);
		ctx.equals("mat2::cols()", mat2(0.0).cols(), 2);

		ctx.equals("mat3::rows()", mat3(0.0).rows(), 3);
		ctx.equals("mat3::cols()", mat3(0.0).cols(), 3);

		ctx.equals("mat4::rows()", mat4(0.0).rows(), 4);
		ctx.equals("mat4::cols()", mat4(0.0).cols(), 4);
	}

	
	// vec.h
	{
		vec4 v = {1.0, -2.0, 3.0, -4.0};
		vec4 w = {1.0, -2.0, 3.0, -4.0};

		vec4 z;
		algebra::vec_zeroes(z);

		ctx.equals("vec_zeroes", linf_norm(z), 0.0);
		ctx.equals("vec_copy", linf_norm(v - w), 0.0);

		ctx.equals("vec::operator+", v + w,  vec4{2.0, -4.0, 6.0, -8.0}, vec4_opt);
		ctx.equals("vec::operator-", v - w, vec4{0.0, 0.0, 0.0, 0.0}, vec4_opt);
		ctx.equals("vec::operator* scalar", v * 2.0, vec4{2.0, -4.0, 6.0, -8.0}, vec4_opt);
		ctx.equals("vec::operator/ scalar", v / 2.0, vec4{0.5, -1.0, 1.5, -2.0}, vec4_opt);

		ctx.equals("vec::dot", v * w, 30.0, 1E-08);
		ctx.equals("vec::norm", v.norm(), std::sqrt(30.0), 1E-08);
		ctx.equals("vec::sqr_norm", v.sqr_norm(), 30.0, 1E-08);
	}


	// mat.h
	{
		mat3 A = {
			{1.0, 2.0, 3.0},
			{4.0, 5.0, 6.0},
			{7.0, 8.0, 9.0}
		};

		mat3 B;
		algebra::mat_zeroes(B);

		ctx.equals("mat_zeroes", linf_norm(B), 0.0);

		mat3 C;
		algebra::mat_copy(C, A);
		ctx.equals("mat_copy", linf_norm(A - C), 0.0);

		algebra::mat_swap_rows(C, 0, 2);
		algebra::mat_swap_cols(C, 0, 2);

		mat3 expected = {
			{9.0, 8.0, 7.0},
			{6.0, 5.0, 4.0},
			{3.0, 2.0, 1.0}
		};

		ctx.equals("mat_swap_rows/mat_swap_cols", linf_norm(C - expected), 0.0);
	}

	{

		vec3 diag = {1.0, 2.0, 3.0};
		mat3 D = algebra::diagonal<mat3>(diag);

		ctx.equals("diagonal", D(0, 0), 1.0);
		ctx.equals("diagonal", D(1, 1), 2.0);
		ctx.equals("diagonal", D(2, 2), 3.0);
	}

	{
		mat2 M = {
			{1.0, 2.0},
			{3.0, 4.0}
		};

		mat2 N2 = {
			{5.0, 6.0},
			{7.0, 8.0}
		};
		
		ctx.equals("mat::operator+", linf_norm((M + N2) - mat2{{6.0, 8.0}, {10.0, 12.0}}), 0.0);
		ctx.equals("mat::operator-", linf_norm((N2 - M) - mat2{{4.0, 4.0}, {4.0, 4.0}}), 0.0);
		ctx.equals("mat::operator* scalar", linf_norm((M * 2.0) - mat2{{2.0, 4.0}, {6.0, 8.0}}), 0.0);
		ctx.equals("mat::transform(vector)", linf_norm((M * vec2{1.0, 2.0}) - vec2{5.0, 11.0}), 0.0);
		ctx.equals("mat::mul(matrix)", linf_norm((M * N2) - mat2{{19.0, 22.0}, {43.0, 50.0}}), 0.0);
	}


	test_residual(ctx, "mat::unpack", [&]() {

		mat<real> A = rand_mat(0.0, 1.0, N, N);
		vec<real> v = A.unpack();
		mat<real> B (v, A.rows(), A.cols());

		return linf_norm(A - B);
	});


	test_residual(ctx, "mat<N, K>::unpack", [&]() {

		mat4 A = rand_mat(0.0, 1.0, 4, 4);
		vec<real, 16> v = A.unpack();
		mat4 B (v, A.rows(), A.cols());

		return linf_norm(A - B);
	});


	// transform.h

	{
		mat3 I3 = algebra::identity<mat3>();
		ctx.equals("identity", linf_norm(I3 - algebra::identity<mat3>()), 0.0);

		vec3 t = {4.0, 5.0, 6.0};
		mat4 T = algebra::translation<mat4>(t);

		ctx.equals("translation", T(0, 3), 4.0);
		ctx.equals("translation", T(1, 3), 5.0);
		ctx.equals("translation", T(2, 3), 6.0);
		ctx.equals("translation", T(3, 3), 1.0);

		const real half_pi = PI / 2.0;
		mat2 R2 = algebra::rotation_2d<mat2>(half_pi);

		ctx.equals("rotation_2d", R2(0, 0), 0.0, 1E-08);
		ctx.equals("rotation_2d", R2(0, 1), -1.0, 1E-08);
		ctx.equals("rotation_2d", R2(1, 0), 1.0, 1E-08);
		ctx.equals("rotation_2d", R2(1, 1), 0.0, 1E-08);

		mat3 Rx = algebra::rotation_3d_xaxis<mat3>(half_pi);
		mat4 Ry = algebra::rotation_3d_yaxis<mat4>(half_pi, 4, 4);
		mat4 Rz = algebra::rotation_3d_zaxis<mat4>(half_pi, 4, 4);

		ctx.equals("rotation_3d_xaxis", Rx(1, 1), 0.0, 1E-08);
		ctx.equals("rotation_3d_xaxis", Rx(1, 2), -1.0, 1E-08);
		ctx.equals("rotation_3d_xaxis", Rx(2, 1), 1.0, 1E-08);
		ctx.equals("rotation_3d_yaxis", Ry(0, 0), 0.0, 1E-08);
		ctx.equals("rotation_3d_yaxis", Ry(0, 2), 1.0, 1E-08);
		ctx.equals("rotation_3d_zaxis", Rz(0, 1), -1.0, 1E-08);
		ctx.equals("rotation_3d_zaxis", Rz(1, 0), 1.0, 1E-08);

		mat4 P = algebra::perspective<mat4>(-1.0, 1.0, -1.0, 1.0, 1.0, 10.0, 4, 4);

		ctx.equals("perspective", P(0, 0), 1.0, 1E-08);
		ctx.equals("perspective", P(1, 1), 1.0, 1E-08);
		ctx.equals("perspective", P(2, 3), -1.0, 1E-08);
	}


	// distance.h

	{
		vec3 a = {1.0, -2.0, 3.0};
		vec3 b = {4.0, 1.0, -1.0};

		ctx.equals("l1_norm", algebra::l1_norm(a), 6.0);
		ctx.equals("l2_norm", algebra::l2_norm(a), std::sqrt(14.0), 1E-08);
		ctx.equals("linf_norm", algebra::linf_norm(a), 3.0);
		ctx.equals("euclidean_distance", algebra::euclidean_distance(a, a), 0.0);
		ctx.equals("manhattan_distance", algebra::manhattan_distance(a, a), 0.0);
		ctx.equals("chebyshev_distance", algebra::chebyshev_distance(a, a), 0.0);
		ctx.equals("minkowski_distance", algebra::minkowski_distance(a, b, 2), algebra::euclidean_distance(a, b), 1E-08);
		ctx.equals("discrete_distance", algebra::discrete_distance(a, a), 0.0);
		ctx.equals("canberra_distance", algebra::canberra_distance(a, a), 0.0);
	}


	// parallel.h

	{
		vec4 v = {1.0, -2.0, 3.0, -4.0};

		vec4 squared = parallel::square(v);
		vec4 rooted = parallel::sqrt(vec4{1.0, 4.0, 9.0, 16.0});
		vec4 mapped = parallel::map([](real x) { return x + 1.0; }, v);
		vec4 transformed = v;

		parallel::transform([](real x) { return x * 2.0; }, transformed);

		ctx.equals("parallel::square", squared,  vec4{1.0, 4.0, 9.0, 16.0}, vec4_opt);
		ctx.equals("parallel::sqrt", rooted, vec4{1.0, 2.0, 3.0, 4.0}, vec4_opt);
		ctx.equals("parallel::map", mapped, vec4{2.0, -1.0, 4.0, -3.0}, vec4_opt);
		ctx.equals("parallel::transform", transformed, vec4{2.0, -4.0, 6.0, -8.0}, vec4_opt);
		ctx.equals("parallel::pow", parallel::pow(vec3{1.0, 2.0, 3.0}, 2), vec3{1.0, 4.0, 9.0}, vec3_opt);
		ctx.equals("parallel::powf", parallel::powf(vec3{1.0, 2.0, 3.0}, 2.0), vec3{1.0, 4.0, 9.0}, vec3_opt);
		ctx.equals("parallel::exp", parallel::exp(vec3{0.0, 1.0, 2.0}), vec3{1.0, std::exp(1.0), std::exp(2.0)}, vec3_opt);
		ctx.equals("parallel::ln", parallel::ln(vec3{1.0, std::exp(1.0), std::exp(2.0)}), vec3{0.0, 1.0, 2.0}, vec3_opt);
	}

}
