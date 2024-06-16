
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"
#include <ctime>

using namespace chebyshev;
using namespace theoretica;

const unsigned int M = 100;


prec::estimate_result test_generator(PRNG g, interval k, Real tol, unsigned int n) {

	Real max = 0;
	Real sum = 0;
	Real sum2 = 0;

	for (size_t j = 0; j < M; ++j) {
		
		real sample_sum = 0;

		for (size_t i = 0; i < n; ++i)
			sample_sum += rand_uniform(k.a, k.b, g);

		real mean = std::abs(sample_sum / n);
		max = th::max(max, mean);
		sum += mean;
		sum += square(mean);
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

	const size_t N = 1000000;

	prec::state.outputFolder = "test/";
	prec::state.defaultIterations = N;

	// Tolerance of 4 standard deviations of the mean
	prec::state.defaultTolerance = (1.0 / SQRT3 / std::sqrt(N)) * 4;
	
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

		// Test autocorrelation of sequences

		vec<> rng_xoshiro(N);
		PRNG g_xoshiro = PRNG::xoshiro(time(nullptr));
		pdf_sampler::uniform(0, 1, g_xoshiro).fill(rng_xoshiro, N);
		prec::equals("PRNG::xoshiro (autocorr.)",
			th::autocorrelation(rng_xoshiro), 0, 1E-2);

		vec<> rng_wyrand(N);
		PRNG g_wyrand = PRNG::wyrand(time(nullptr));
		pdf_sampler::uniform(0, 1, g_wyrand).fill(rng_wyrand, N);
		prec::equals("PRNG::wyrand (autocorr.)",
			th::autocorrelation(rng_wyrand), 0, 1E-2);

		vec<> rng_lincong(N);
		PRNG g_lincong = PRNG::linear_congruential(time(nullptr));
		pdf_sampler::uniform(0, 1, g_lincong).fill(rng_lincong, N);
		prec::equals("PRNG::linear_congruential (autocorr.)",
			th::autocorrelation(rng_lincong), 0, 1E-2);

		vec<> rng_splitmix(N);
		PRNG g_splitmix = PRNG::splitmix64(time(nullptr));
		pdf_sampler::uniform(0, 1, g_splitmix).fill(rng_splitmix, N);
		prec::equals("PRNG::splitmix64 (autocorr.)",
			th::autocorrelation(rng_splitmix), 0, 1E-2);

		vec<> rng_middlesquare(N);
		PRNG g_middlesquare = PRNG::middlesquare(time(nullptr));
		pdf_sampler::uniform(0, 1, g_middlesquare).fill(rng_middlesquare, N);
		prec::equals("PRNG::middlesquare (autocorr.)",
			th::autocorrelation(rng_middlesquare), 0, 1E-2);


		// Metropolis
		{
			PRNG g = PRNG::xoshiro(time(nullptr));
			pdf_sampler gauss = pdf_sampler::gaussian(0, 1, g);

			std::vector<real> sample;
			real res = 0;

			for (int i = 0; i < 100; ++i) {
				for (int j = 0; j < 1000; ++j) {
					sample.push_back(
						metropolis(
							[](real x) {
								
								if(x < 0)
									return 0.0;

								return th::exp(-x);
							}, gauss, 1));
				}

				res += mean(sample);
				sample.clear();
			}

			res /= 100;
			prec::equals("metropolis", res, 1, 0.05);
		}

	prec::terminate();
}
