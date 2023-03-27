
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev/prec.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


// Generate a random matrix with elements inside the interval [a, b]
template<unsigned int N, unsigned int M>
mat<N, M> rand_mat(real a, real b, PRNG& g) {

	mat<N, M> A;
	for (size_t i = 0; i < N; ++i)
		for (size_t j = 0; j < M; ++j)
			A.iat(i, j) = rand_uniform(a, b, g);

	return A;
}


// Test mat<N, N>::inverse()
template<unsigned int N>
prec::estimate_result test_matrix_inverse(interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (size_t i = 0; i < n; ++i) {

		mat<N, N> A = rand_mat<N, N>(k.a, k.b, g);

		// Skip singular matrices
		if(th::abs(A.det()) <= MACH_EPSILON) {
			if(i) i--;
			continue;
		}

		// Resulting matrix expected to be identity
		mat<N, N> R = A * A.inverse();

		for (size_t j = 0; j < N; ++j) {
			for (size_t k = 0; k < N; ++k) {

				Real diff = th::abs(R.iat(j, k) - kronecker_delta(j, k));

				sum += diff;
				sum2 += square(diff);
				if(diff > max)
					max = diff;
			}
		}

	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum / n;
	res.rms_err = th::sqrt(sum2) / n;
	res.mean_err = sum / N / n;

	// Undefined relative error
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}


// Test mat<N, N>::det()
template<unsigned int N>
prec::estimate_result test_matrix_det(interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (size_t i = 0; i < n; ++i) {

		mat<N, N> A = mat<N, N>();
		for (size_t j = 0; j < N; ++j)
			A.at(j, j) = rand_uniform(k.a, k.b, g);

		real expected = 1;
		for (size_t j = 0; j < N; ++j)
			expected *= A.at(j, j);
		
		real computed = A.det();
		real diff = th::abs(computed - expected);
		
		sum += diff;
		sum2 += square(diff);
		if(diff > max)
			max = diff;
	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum / n;
	res.rms_err = th::sqrt(sum2) / n;
	res.mean_err = sum / N / n;

	// Undefined relative error
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}


prec::estimate_result test_matrix_mul(interval k, Real tol, unsigned int n) {


	mat<3, 3> A = {
		{1, 5, 9},
		{10, 7, 18},
		{3, 11, 5}
	};


	mat<3, 4> B = {
		{7, 5, 0, 11},
		{4, 12, 1, 6},
		{3, 7, 9, 0},
	};

	mat<3, 4> C = {
		{54, 128, 86, 41},
		{152, 260, 169, 152},
		{80, 182, 56, 99}
	};


	mat<3, 4> res = A * B;

	real max = 0;
	real sum = 0;
	real sum2 = 0;

	for (unsigned int i = 0; i < res.row_size(); ++i) {
		for (unsigned int j = 0; j < res.col_size(); ++j) {
			
			real diff = th::abs(res.iat(i, j) - C.iat(i, j));

			if(max < diff)
				max = diff;

			sum += diff;
			sum2 += square(diff);
		}
	}

	prec::estimate_result p;
	p.max_err = max;
	p.abs_err = sum / res.size();
	p.rms_err = th::sqrt(sum2) / res.size();
	p.mean_err = sum / res.size();

	// Undefined relative error
	p.rel_err = 0;

	if(p.max_err > tol)
		p.failed = true;

	return p;
}


template<unsigned int N>
prec::estimate_result test_distance(
	real(*d)(vec<N, real>, vec<N, real>), interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (size_t i = 0; i < n; ++i) {

		vec<N, real> v;

		for (unsigned int j = 0; j < N; ++j) {
			v[j] = rand_uniform(k.a, k.b, g);
		}

		vec<N, real> w = v;		
		real diff = th::abs(d(v, w));
		
		sum += diff;
		sum2 += square(diff);
		if(diff > max)
			max = diff;
	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum / n;
	res.rms_err = th::sqrt(sum2) / n;
	res.mean_err = sum / N / n;

	// Undefined relative error
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";

	std::vector<interval> intervals = {
		interval(-1, 1),
		interval(-1000000, 1000000)
	};

	prec::setup("algebra", argc, argv);

		// Matrices (mat.h)

		prec::estimate("mat2::inverse", test_matrix_inverse<2>, intervals);
		prec::estimate("mat3::inverse", test_matrix_inverse<3>, intervals);
		prec::estimate("mat4::inverse", test_matrix_inverse<4>, intervals);
		prec::estimate("mat10::inverse", test_matrix_inverse<10>, intervals);

		prec::estimate("mat2::det", test_matrix_det<2>, intervals);
		prec::estimate("mat3::det", test_matrix_det<3>, intervals);
		prec::estimate("mat4::det", test_matrix_det<4>, intervals);
		prec::estimate("mat10::det", test_matrix_det<10>, intervals);

		prec::estimate("mat3::operator*", test_matrix_mul, interval(0, 1));

		// Distances and norms (distance.h)

		// Test Lp norms from 1 to 10
		for (unsigned int p = 1; p <= 10; ++p)
			prec::equals("lp_norm<vec3>", lp_norm(vec<3>(0), p), 0);

		prec::equals("lp_norm<vec100>", lp_norm(vec<100>(0), 2), 0);

		// L1
		prec::equals("l1_norm<vec3>", l1_norm(vec<3>(0)), 0);
		prec::equals("l1_norm<vec100>", l1_norm(vec<100>(0)), 0);

		prec::equals("l1_norm<vec4>", l1_norm(vec<4>(1)), 4);
		prec::equals("l1_norm<vec100>", l1_norm(vec<100>(1)), 100);

		// L2
		prec::equals("l2_norm<vec3>", l2_norm(vec<3>(0)), 0);
		prec::equals("l2_norm<vec100>", l2_norm(vec<100>(0)), 0);

		prec::equals("l2_norm<vec4>", l2_norm(vec<4>(1)), 2);
		prec::equals("l2_norm<vec9>", l2_norm(vec<9>(1)), 3);

		// Linf
		prec::equals("linf_norm<vec3>", linf_norm(vec<3>(0)), 0);
		prec::equals("linf_norm<vec100>", linf_norm(vec<100>(0)), 0);
		prec::equals("linf_norm<vec100>", linf_norm(vec<100>(1)), 1);

		// Distances
		prec::estimate("euclidean_distance<3>",
			[](interval k, Real tol, unsigned int n) {
				return test_distance<3>(euclidean_distance<vec<3>>, k, tol, n);
			}, intervals);


		prec::estimate("manhattan_distance<3>",
			[](interval k, Real tol, unsigned int n) {
				return test_distance<3>(manhattan_distance<vec<3>>, k, tol, n);
			}, intervals);

	prec::terminate();
}
