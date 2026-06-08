
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


real f(vec<real, 2> x) {
	return th::sqrt(x * x);
}

real g(real x) {
    return th::sqrt(x);
}


int main(int argc, char const *argv[]) {
	
	auto ctx = benchmark::make_context("montecarlo", argc, argv);
	ctx.output->settings.outputFiles = { "test/benchmark/benchmark_montecarlo.csv" };
	ctx.settings.defaultRuns = 4;
	ctx.settings.multithreading = false;

	th::random::XoshiroPrng g (time(nullptr));

	vec<vec2> domain (10, vec2({0.0, 1.0}));
	vec<size_t> Npoints = {1'000, 10'000, 100'000, 1'000'000, 10'000'000};


	for	(auto sz : Npoints) {

		ctx.benchmark(
			"integral_mc (" + std::to_string(sz) + " points)",
			[&](real x) {
				return integral_mc(f, domain, g, sz).value;
			},
			std::vector<real>(1)
		);
	}

}