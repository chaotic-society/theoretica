
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


// Compare to standard for loop
template<typename Function>
auto wrap_base(Function f) {

	return [f](const vec<real>& v) {

		vec<real> x = v;

		for (real& val : x)
			val = f(val);

		return x[0];
	};
}


int main(int argc, char const *argv[]) {


	// Vector size
	const size_t N = 1'000'000;
	
	// Number of vectors per run
	const size_t M = 10;

	
	auto ctx = benchmark::make_context("parallel", argc, argv);

	ctx.output->settings.outputFiles = { "test/benchmark/benchmark_parallel.csv" };
	ctx.settings.defaultRuns = 10;
	
	th::random::XoshiroPrng g (time(nullptr));
	auto unif = th::random::PdfSampler<th::random::XoshiroPrng>::uniform(0.0, 10.0, g);

	// Generate a uniform sample
	std::vector<vec<real>> data (M, vec<real>(N));

	for (auto& v : data)
		unif.fill(v);

	// To be compared to benchmark_real_analysis,
	// with time multiplied by N.
	// For a function taking t = 1 x 1E-05 ms,
	// N * t = 100 ms

	ctx.benchmark(
		"th::sqrt",
		wrap_base(th::sqrt<real>),
		data
	);

	ctx.benchmark(
		"parallel::sqrt",
		[](const vec<real>& v){ return th::sqrt(v)[0]; },
		data
	);

	ctx.benchmark(
		"th::exp",
		wrap_base(th::exp<real>),
		data
	);
	
	ctx.benchmark(
		"parallel::exp",
		[](const vec<real>& v){ return th::exp(v)[0]; },
		data
	);

	ctx.benchmark(
		"th::atan",
		wrap_base(th::atan<real>),
		data
	);

	ctx.benchmark(
		"parallel::atan", 
		[](const vec<real>& v){ return th::atan(v)[0]; },
		data
	);
}
