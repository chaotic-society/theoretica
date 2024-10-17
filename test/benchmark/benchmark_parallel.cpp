
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;


template<typename Function>
auto wrap(Function f) {

	return [f](const vec<real>& v) {
		return f(v)[0];
	};
}


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

	benchmark::setup("parallel", argc, argv);

		// Vector size
		const size_t N = 1'000'000;
		
		// Number of vectors per run
		const size_t M = 10;

		output::settings.outputFiles = { "test/benchmark/benchmark_parallel.csv" };
		benchmark::settings.defaultRuns = 10;

		PRNG g = PRNG::xoshiro(time(nullptr));
		pdf_sampler unif = pdf_sampler::uniform(0.0, 10.0, g);

		// Generate a uniform sample
		std::vector<vec<real>> data (M, vec<real>(N));

		for (auto& v : data)
			unif.fill(v);

		// To be compared to benchmark_real_analysis,
		// with time multiplied by N.
		// For a function taking t = 1 x 1E-05 ms,
		// N * t = 100 ms

		benchmark::benchmark(
			"base_sqrt",
			wrap_base(th::sqrt<real>),
			data
		);

		benchmark::benchmark(
			"parallel::sqrt",
			wrap(parallel::sqrt<vec<real>>),
			data
		);

		benchmark::benchmark(
			"parallel::square",
			wrap(parallel::square<vec<real>>),
			data
		);

		benchmark::benchmark(
			"base_exp",
			wrap_base(th::exp<real>),
			data
		);
		
		benchmark::benchmark(
			"parallel::exp",
			wrap(parallel::exp<vec<real>>),
			data
		);

		benchmark::benchmark(
			"base_atan",
			wrap_base(th::atan<real>),
			data
		);

		benchmark::benchmark(
			"parallel::atan", 
			wrap(parallel::atan<vec<real>>),
			data
		);

	benchmark::terminate();
}
