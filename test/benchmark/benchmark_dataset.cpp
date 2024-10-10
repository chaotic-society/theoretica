
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	benchmark::setup("dataset", argc, argv);

		output::state.outputFiles = { "test/benchmark/benchmark_dataset.csv" };

		const size_t N = 1'000'000;
		auto dummy = [](unsigned int i) { return 0; };

		auto opt = benchmark::benchmark_options<real>(
			10, 10, dummy
		);

		PRNG g = PRNG::xoshiro(time(nullptr));
	    pdf_sampler gauss = pdf_sampler::gaussian(0, 1'000'000, g);

	    // Generate a gaussian sample
	    vec<real> v = vec<real>(N);
	    gauss.fill(v);

	    benchmark::benchmark(
	    	"sum",
	    	[v](real x) { return sum(v); },
	    	opt
	    );

		benchmark::benchmark(
			"sum_pairwise",
			[v](real x) { return sum_pairwise(v); },
			opt
		);

		benchmark::benchmark(
			"sum_compensated",
			[v](real x) { return sum_compensated(v); },
			opt
		);

	benchmark::terminate();
}
