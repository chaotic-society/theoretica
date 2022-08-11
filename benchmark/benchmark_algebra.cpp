
#include "theoretica.h"
#include "chebyshev/benchmark.h"
using namespace chebyshev;
using namespace theoretica;
using benchmark::benchmark_result;

#include <ctime>


template<unsigned int N>
benchmark_result benchmark_mat_inverse(unsigned int iter, unsigned int runs) {

	std::vector<mat<N, N>> A;
	A.reserve(iter);
	real c = 0;

	const unsigned int MAX = 100000;
	const unsigned int MIN = -100000;

	PRNG g = PRNG::xoshiro(time(nullptr));

	// Generate random invertible matrices
	for (size_t i = 0; i < iter; ++i) {
		do {
			for (size_t j = 0; j < N; ++j) {
				for (size_t k = 0; k < N; ++k) {
					A[i].at(j, k) = rand_uniform(MIN, MAX, g);
				}
			}
		} while(A[i].det() < MACH_EPSILON);
	}

	long double elapsed = 0;

	for (size_t i = 0; i < iter; ++i) {
		
		timer t = timer();

		for (size_t i = 0; i < iter; ++i)
			c += A[i].inverse().get(0, 0);

		elapsed += t();
	}

	return benchmark_result(elapsed, iter, runs);
}



int main(int argc, char const *argv[]) {

	benchmark::state.outputFolder = "benchmark/";
	
	benchmark::setup("algebra", argc, argv);

		benchmark::custom_request("mat2::inverse()", benchmark_mat_inverse<2>, 1000, 1000);
		benchmark::custom_request("mat3::inverse()", benchmark_mat_inverse<3>, 1000, 1000);
		benchmark::custom_request("mat4::inverse()", benchmark_mat_inverse<4>, 1000, 1000);
		benchmark::custom_request("mat10::inverse()", benchmark_mat_inverse<10>, 1000, 1000);

	benchmark::terminate();
}
