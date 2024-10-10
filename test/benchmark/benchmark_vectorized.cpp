
#include "theoretica.h"
#include "chebyshev/benchmark.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	const size_t N = 1'000'000;

	auto dummy = [](unsigned int i) { return 0; };

	benchmark::state.outputFolder = "test/benchmark/";
	benchmark::state.defaultIterations = 100;
	benchmark::state.defaultRuns = 1;


	benchmark::setup("vectorized", argc, argv);


		PRNG g = PRNG::xoshiro(time(nullptr));
		pdf_sampler unif = pdf_sampler::uniform(0, 10, g);

		// Generate a uniform sample
		vec<real> v = vec<real>(N);
		unif.fill(v);


		// To be compared to benchmark_real_analysis,
		// with time multiplied by N.
		// For a function taking t = 1 x 1E-05 ms,
		// N * t = 100 ms

		benchmark::request("vectorized::square",
			[v](real x) { return vectorized::square(v)[0]; }, dummy);

		benchmark::request("vectorized::sqrt",
			[v](real x) { return vectorized::sqrt(v)[0]; }, dummy);

		benchmark::request("vectorized::exp",
			[v](real x) { return vectorized::exp(v)[0]; }, dummy);

		benchmark::request("vectorized::pow(10)",
			[v](real x) { return vectorized::pow(v, 10)[0]; }, dummy);

		benchmark::request("vectorized::atan",
			[v](real x) { return vectorized::atan(v)[0]; }, dummy);

	benchmark::terminate();
}
