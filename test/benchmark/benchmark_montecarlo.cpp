
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


real f(vec<real> x) {
	return th::sqrt(x * x);
}


int main(int argc, char const *argv[]) {
	
	auto ctx = benchmark::make_context("montecarlo", argc, argv);
	ctx.output->settings.outputFiles = { "test/benchmark/benchmark_montecarlo.csv" };
	ctx.settings.defaultRuns = 10;
	ctx.settings.multithreading = false;

	PRNG g1 = PRNG::xoshiro(13631615);
	PRNG g2 = PRNG::xoshiro(16473473);

	vec<vec2> domain (10, vec2({0.0, 1.0}));
	vec<size_t> Npoints = {100'000, 1'000'000, 10'000'000};


	for	(auto sz : Npoints) {

		ctx.benchmark(
			"integral_crude (" + std::to_string(sz) + " points)",
			[&](real x) {
				return integral_crude(f, domain, g1, sz);
			},
			std::vector<real>(1), 1
		);

		ctx.benchmark(
			"integral_mc (" + std::to_string(sz) + " points)",
			[&](real x) {
				return integral_mc(f, domain, g2, sz);
			},
			std::vector<real>(1), 1
		);
	}

}
