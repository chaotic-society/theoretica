
#include "theoretica.h"
#include "chebyshev.h"
using namespace chebyshev;
using namespace theoretica;

// A macro to benchmark a real function
#define BENCHMARK_REAL(func, opt) \
	benchmark::benchmark(#func, CAST_LAMBDA(func, real), opt)


int main(int argc, char const *argv[]) {
	
	benchmark::setup("real_analysis", argc, argv);

		output::state.outputFolder = "benchmark/";

		// Benchmark options for real functions
		auto R_opt = benchmark::benchmark_options<real>(
			10, 100'000,
			benchmark::generator::uniform1D(-1E+06, +1E+06)
		);

		// Benchmark options for functions over the positive reals
		auto Rplus_opt = benchmark::benchmark_options<real>(
			10, 100'000,
			benchmark::generator::uniform1D(0, +1E+06)
		);

		BENCHMARK_REAL(th::sqrt, R_opt);

		BENCHMARK_REAL(th::square, R_opt);
		BENCHMARK_REAL(th::cube, R_opt);
		BENCHMARK_REAL(th::isqrt<uint32_t>, Rplus_opt);
		BENCHMARK_REAL(th::icbrt<uint32_t>, Rplus_opt);
		BENCHMARK_REAL(th::sqrt, Rplus_opt);
		BENCHMARK_REAL(th::cbrt, R_opt);
		BENCHMARK_REAL(th::abs, R_opt);
		BENCHMARK_REAL(th::sgn, R_opt);
		BENCHMARK_REAL(th::floor, R_opt);
		BENCHMARK_REAL(th::fract, R_opt);

		benchmark::benchmark("th::max",
			[](real x) { return max(-1E-09, x); },
			R_opt
		);

		benchmark::benchmark("th::max",
			[](real x) { return max(x, 1E+09); },
			R_opt
		);

		benchmark::benchmark("th::min",
			[](real x) { return min(-1E-09, x); },
			R_opt
		);

		benchmark::benchmark("th::min",
			[](real x) { return min(x, 1E+09); },
			R_opt
		);

		benchmark::benchmark("th::clamp",
			[](real x) { return clamp(x, -1E+09, +1E+09); },
			R_opt
		);

		benchmark::benchmark("th::clamp",
			[](real x) { return clamp(x, 0, 1); },
			R_opt
		);

		BENCHMARK_REAL(th::ln, Rplus_opt);
		BENCHMARK_REAL(th::log2, Rplus_opt);
		BENCHMARK_REAL(th::log10, Rplus_opt);

		benchmark::benchmark<real>(
			"th::pow",
			[](real x) { return th::pow(1.1, (int) x); },
			10, 100'000, benchmark::generator::uniform1D(-100.0, 100.0)
		);

		benchmark::benchmark<real>(
			"th::pow",
			[](real x) { return th::pow(1.1, (int) -x); },
			10, 100'000, benchmark::generator::uniform1D(-100.0, 100.0)
		);

		benchmark::benchmark(
			"th::root(x, 10)",
			[](real x) { return th::root(x, 10); },
			Rplus_opt
		);

		BENCHMARK_REAL(th::fract, R_opt);
		
		BENCHMARK_REAL(
			th::exp,
			benchmark::benchmark_options<real>(
				10, 100'000, benchmark::generator::uniform1D(-100, 10)
			)
		);

		benchmark::benchmark(
			"th::powf (1)",
			[](real x) { return th::powf(x, 10); },
			R_opt
		);

		benchmark::benchmark(
			"th::powf (2)",
			[](real x) { return th::powf(x, -10); },
			R_opt
		);

		BENCHMARK_REAL(th::sin, R_opt);
		BENCHMARK_REAL(th::cos, R_opt);
		BENCHMARK_REAL(th::tan, R_opt);
		BENCHMARK_REAL(th::cot, R_opt);

		BENCHMARK_REAL(th::atan, R_opt);
		BENCHMARK_REAL(th::asin, R_opt);
		BENCHMARK_REAL(th::acos, R_opt);

		auto exp_opt = benchmark::benchmark_options<real>(
			10, 100'000,
			benchmark::generator::uniform1D(-50, 50)
		);

		BENCHMARK_REAL(th::sinh, exp_opt);
		BENCHMARK_REAL(th::cosh, exp_opt);
		BENCHMARK_REAL(th::tanh, exp_opt);
		BENCHMARK_REAL(th::coth, exp_opt);
		BENCHMARK_REAL(th::sigmoid, exp_opt);
		BENCHMARK_REAL(th::sinc, R_opt);
		BENCHMARK_REAL(th::heaviside, R_opt);
		BENCHMARK_REAL(th::radians, R_opt);
		BENCHMARK_REAL(th::degrees, R_opt);

		benchmark::benchmark(
			"th::binomial_coeff",
			[](real x) { return binomial_coeff<uint64_t>(20, x); },
			benchmark::benchmark_options<real>(
				10, 10'000, benchmark::generator::uniform1D(0, 20)
			)
		);

		BENCHMARK_REAL(
			th::special::gamma,
			benchmark::benchmark_options<real>(
				10, 10'000, benchmark::generator::uniform1D(0.1, 10)
			)
		);

	benchmark::terminate();
}
