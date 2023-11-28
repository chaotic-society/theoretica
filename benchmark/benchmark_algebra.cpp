
#include "theoretica.h"
#include "chebyshev/benchmark.h"
#include <ctime>

using namespace chebyshev;
using namespace theoretica;
using benchmark::benchmark_result;


template<unsigned int N, unsigned int M = N>
mat<real, N, M> rand_mat(real min, real max, PRNG& g) {

	mat<real, N, M> A;

	for (unsigned int i = 0; i < N; ++i)
		for (unsigned int j = 0; j < M; ++j)
			A.at(i, j) = rand_uniform(min, max, g);

	return A;
}


template<unsigned int N>
benchmark_result benchmark_mat_det(unsigned int iter, unsigned int runs) {

	std::vector<mat<real, N, N>> A(iter);
	volatile real c = 0;

	const int MAX = 100000;
	const int MIN = -100000;
	PRNG g = PRNG::xoshiro(time(nullptr));

	for (unsigned int i = 0; i < iter; ++i)
		A[i] = rand_mat<N>(MIN, MAX, g);

	long double elapsed = 0;

	for (unsigned int i = 0; i < runs; ++i) {
		
		timer t = timer();

		for (unsigned int j = 0; j < iter; ++j)
			c += A[j].det();

		elapsed += t();
	}

	return benchmark_result(elapsed, iter, runs);
}


template<unsigned int N>
benchmark_result benchmark_mat_inverse(unsigned int iter, unsigned int runs) {

	std::vector<mat<real, N, N>> A(iter);
	volatile real c = 0;

	const int MAX = 100000;
	const int MIN = -100000;
	PRNG g = PRNG::xoshiro(time(nullptr));

	// Generate random invertible matrices
	for (unsigned int i = 0; i < iter; ++i) {
		do {
			A[i] = rand_mat<N>(MIN, MAX, g);
		} while(A[i].det() < MACH_EPSILON);
	}

	long double elapsed = 0;

	for (unsigned int i = 0; i < runs; ++i) {
		
		timer t = timer();

		for (unsigned int j = 0; j < iter; ++j)
			c += A[j].inverse().get(0, 0);

		elapsed += t();
	}

	return benchmark_result(elapsed, iter, runs);
}



int main(int argc, char const *argv[]) {

	benchmark::state.outputFolder = "benchmark/";
	
	benchmark::setup("algebra", argc, argv, 1000, 1000);

		benchmark::custom_request("mat2::inverse()", benchmark_mat_inverse<2>);
		benchmark::custom_request("mat3::inverse()", benchmark_mat_inverse<3>);
		benchmark::custom_request("mat4::inverse()", benchmark_mat_inverse<4>);
		benchmark::custom_request("mat10::inverse()", benchmark_mat_inverse<10>);

		benchmark::custom_request("mat2::det()", benchmark_mat_det<2>);
		benchmark::custom_request("mat3::det()", benchmark_mat_det<3>);
		benchmark::custom_request("mat4::det()", benchmark_mat_det<4>);
		benchmark::custom_request("mat10::det()", benchmark_mat_det<10>);

	benchmark::terminate();
}
