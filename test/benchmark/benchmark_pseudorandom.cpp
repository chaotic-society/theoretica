
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	// const size_t N = 1000000;
	// const uint64_t seed = time(nullptr);
	// auto dummy = [](unsigned int i) { return 0; };

	benchmark::settings.outputFiles = { "test/benchmark/" };
	benchmark::settings.defaultIterations = 10;
	benchmark::settings.defaultRuns = 10;

	benchmark::setup("pseudorandom", argc, argv);

		// PRNG g_xoshiro = PRNG::xoshiro(seed);
		// PRNG g_wyrand = PRNG::wyrand(seed);
		// PRNG g_lc = PRNG::linear_congruential(seed);
		// PRNG g_splitmix64 = PRNG::splitmix64(seed);
		// PRNG g_middlesquare = PRNG::middlesquare(seed);

		// // Allocate a vector to fill with random numbers
		// vec<uint64_t> v = vec<uint64_t>(N);

		// Measure the time taken to generate 1 million numbers

		// benchmark::request("PRNG::xoshiro",
		// 	[&](real x) {
		// 		for (size_t i = 0; i < N; ++i)
		// 			v[i] = g_xoshiro();
		// 		return v[0];
		// 	}, dummy);

		// benchmark::request("PRNG::wyrand",
		// 	[&](real x) {
		// 		for (size_t i = 0; i < N; ++i)
		// 			v[i] = g_wyrand();
		// 		return v[0];
		// 	}, dummy);

		// benchmark::request("PRNG::linear_...",
		// 	[&](real x) {
		// 		for (size_t i = 0; i < N; ++i)
		// 			v[i] = g_lc();
		// 		return v[0];
		// 	}, dummy);

		// benchmark::request("PRNG::splitmix64",
		// 	[&](real x) {
		// 		for (size_t i = 0; i < N; ++i)
		// 			v[i] = g_splitmix64();
		// 		return v[0];
		// 	}, dummy);

		// benchmark::request("PRNG::middlesquare",
		// 	[&](real x) {
		// 		for (size_t i = 0; i < N; ++i)
		// 			v[i] = g_middlesquare();
		// 		return v[0];
		// 	}, dummy);

	benchmark::terminate();
}
