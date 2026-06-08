
#include <ctime>

#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


// Sample size for PRNGs
const size_t N = 1'000'000;

// Iterations
const unsigned int iterations = 1000;

// Runs
const unsigned int runs = 1;


// Benchmark a pseudorandom number generator
template<typename Prng>
void benchmark_PRNG(std::string PRNGName, benchmark::benchmark_context& ctx, Prng& g) {

	auto opt = benchmark::benchmark_options<real>(iterations, runs);

	ctx.benchmark(PRNGName + " (1M)",
		[&](uint64_t x) {

			for (size_t i = 0; i < N - 1; ++i)
				g();

			return g();
	}, opt);
}


int main(int argc, char const *argv[]) {

	auto ctx = benchmark::make_context("random", argc, argv);

	ctx.settings.outputFiles = { "test/benchmark/benchmark_random.csv" };

	// Disable multithreading to avoid segmentation faults caused by PRNG benchmark
	ctx.settings.multithreading = false;

	const uint64_t seed = time(nullptr);
	th::random::SplitmixPrng splitmix (seed);
	th::random::XoshiroPrng xoshiro (seed);
	th::random::WyrandPrng wyrand (seed);

	// Measure the time taken to generate 1 million numbers
	benchmark_PRNG("random::XoshiroPrng", ctx, xoshiro);
	benchmark_PRNG("random::WyrandPrng", ctx, wyrand);
	benchmark_PRNG("random::SplitmixPrng", ctx, splitmix);
}
