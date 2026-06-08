
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	auto ctx = benchmark::make_context("dataset", argc, argv);
	ctx.output->settings.outputFiles = { "test/benchmark/benchmark_dataset.csv" };

	// Disable multithreading because it would
	// interfere with benchmarks
	ctx.settings.multithreading = false;

	const size_t N = 1'000'000;
	auto opt = benchmark::benchmark_options<real>(10, 10);

	th::random::XoshiroPrng g (time(nullptr));
	auto gauss = th::random::PdfSampler<th::random::XoshiroPrng>::gaussian(0, 1'000'000, g);

	// Generate a gaussian sample
	vec<real> v = vec<real>(N);
	gauss.fill(v);

	ctx.benchmark(
		"sum",
		[v](real x) { return sum(v); },
		opt
	);

	ctx.benchmark(
		"sum_pairwise",
		[v](real x) { return sum_pairwise(v); },
		opt
	);

	ctx.benchmark(
		"sum_compensated",
		[v](real x) { return sum_compensated(v); },
		opt
	);
}
