
#include "theoretica.h"
#include "chebyshev/benchmark.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	const size_t N = 1000000;

	auto dummy = [](unsigned int i) { return 0; };

	benchmark::state.outputFolder = "benchmark/";
	benchmark::state.defaultIterations = 10;
	benchmark::state.defaultRuns = 10;

	benchmark::setup("real_analysis", argc, argv);


		PRNG g = PRNG::xoshiro(time(nullptr));
	    pdf_sampler gauss = pdf_sampler::gaussian(0, 1000000, g);

	    // Generate a gaussian sample
	    vec<real> v = vec<real>(N);
	    gauss.fill(v, N);

	    benchmark::request("sum",
			[v](real x) { return sum(v); }, dummy);

		benchmark::request("sum_pairwise",
			[v](real x) { return sum_pairwise(v); }, dummy);

		benchmark::request("sum_compensated",
			[v](real x) { return sum_compensated(v); }, dummy);

	benchmark::terminate();
}
