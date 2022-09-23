
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



int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";

	std::vector<interval> intervals = {
		interval(-1, 1),
		interval(-1000000, 1000000)
	};

	prec::setup("algebra", argc, argv);

		prec::estimate("mat2::inverse", test_matrix_inverse<2>, intervals);
		prec::estimate("mat3::inverse", test_matrix_inverse<3>, intervals);
		prec::estimate("mat4::inverse", test_matrix_inverse<4>, intervals);
		prec::estimate("mat10::inverse", test_matrix_inverse<10>, intervals);

		prec::estimate("mat2::det", test_matrix_det<2>, intervals);
		prec::estimate("mat3::det", test_matrix_det<3>, intervals);
		prec::estimate("mat4::det", test_matrix_det<4>, intervals);
		prec::estimate("mat10::det", test_matrix_det<10>, intervals);

		prec::estimate("mat3::operator*", test_matrix_mul, interval(0, 1));

	prec::terminate();
}
