
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"
#include <ctime>

using namespace chebyshev;
using namespace theoretica;

const unsigned int M = 10;


prec::estimate_result test_generator(PRNG g, interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	for (size_t j = 0; j < M; ++j) {
		
		real sample_sum = 0;

		for (size_t i = 0; i < n; ++i)
			sample_sum += rand_uniform(k.a, k.b, g);

		real diff = th::abs((sample_sum / n));
		max = th::max(max, diff);
		sum += diff;
		sum += square(diff);
	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum;
	res.rms_err = th::sqrt(sum2) / M;
	res.mean_err = sum / M;

	// Undefined relative error
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}



int main(int argc, char const *argv[]) {

	// Normal interval
	interval I = interval(-1, 1);

	prec::state.outputFolder = "test/";
	prec::state.defaultIterations = 1000000;
	prec::state.defaultTolerance = 0.005;
	
	prec::setup("pseudorandom");

		prec::estimate("PRNG::xoshiro",
			[](interval k, Real tol, unsigned int n) {
				return test_generator(PRNG::xoshiro(time(nullptr)), k, tol, n);
			}, I);

		prec::estimate("PRNG::wyrand",
			[](interval k, Real tol, unsigned int n) {
				return test_generator(PRNG::wyrand(time(nullptr)), k, tol, n);
			}, I);

		prec::estimate("PRNG::linear_congruential",
			[](interval k, Real tol, unsigned int n) {
				return test_generator(PRNG::linear_congruential(time(nullptr)), k, tol, n);
			}, I);

		prec::estimate("PRNG::splitmix64",
			[](interval k, Real tol, unsigned int n) {
				return test_generator(PRNG::splitmix64(time(nullptr)), k, tol, n);
			}, I);

		prec::estimate("PRNG::middlesquare",
			[](interval k, Real tol, unsigned int n) {
				return test_generator(PRNG::middlesquare(time(nullptr)), k, tol, n);
			}, I);

	prec::terminate();
}
