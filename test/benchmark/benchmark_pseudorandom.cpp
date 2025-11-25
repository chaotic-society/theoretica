
#include <ctime>

#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


// Sample size for PRNGs
const size_t N = 1'000'000;


// Benchmark a pseudorandom number generator
void benchmark_PRNG(std::string PRNGName, benchmark::benchmark_context& ctx, PRNG& g) {

	vec<uint64_t> v (N);
	auto opt = benchmark::benchmark_options<real>(50, 1);

	ctx.benchmark("PRNG::" + PRNGName + " (1M)",
		[&](real x) {
			for (size_t i = 0; i < N; ++i)
				v[i] = g();
			return v[0];
	}, opt);
}


int main(int argc, char const *argv[]) {

	auto ctx = benchmark::make_context("pseudorandom", argc, argv);

	ctx.settings.outputFiles = { "test/benchmark/benchmark_pseudorandom.csv" };
	ctx.settings.defaultIterations = 10;
	ctx.settings.defaultRuns = 10;

	// Disable multithreading to avoid segmentation faults
	// caused by PRNG benchmark
	ctx.settings.multithreading = false;

	const uint64_t seed = time(nullptr);
	PRNG g_xoshiro = PRNG::xoshiro(seed);
	PRNG g_wyrand = PRNG::wyrand(seed);
	PRNG g_lc = PRNG::linear_congruential(seed);
	PRNG g_splitmix64 = PRNG::splitmix64(seed);
	PRNG g_middlesquare = PRNG::middlesquare(seed);

	// Allocate a vector to fill with random numbers
	vec<uint64_t> v = vec<uint64_t>(N);

	// Measure the time taken to generate 1 million numbers

	benchmark_PRNG("xoshiro", ctx, g_xoshiro);
	benchmark_PRNG("wyrand", ctx, g_wyrand);
	benchmark_PRNG("linear_congruential", ctx, g_lc);
	benchmark_PRNG("splitmix64", ctx, g_splitmix64);
	benchmark_PRNG("middlesquare", ctx, g_middlesquare);
}
