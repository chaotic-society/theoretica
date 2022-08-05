
#include "theoretica.h"
#include "chebyshev/benchmark.h"
using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	const real MIN = -1000000;
	const real MAX = 1000000;

	benchmark::state.outputFolder = "benchmark/";
	
	benchmark::setup("real_analysis");

		BENCHMARK(th::square, MIN, MAX);
		BENCHMARK(th::cube, MIN, MAX);
		BENCHMARK(th::isqrt<uint32_t>, 0, MAX);
		BENCHMARK(th::icbrt<uint32_t>, 0, MAX);
		BENCHMARK(th::sqrt, 0, MAX);
		BENCHMARK(th::cbrt, MIN, MAX);
		BENCHMARK(th::abs, MIN, MAX);
		BENCHMARK(th::sgn, MIN, MAX);
		BENCHMARK(th::floor, MIN, MAX);
		BENCHMARK(th::fract, MIN, MAX);

		benchmark::request("th::max (1)",
			[MIN, MAX](real x) { return max(MIN, x); }, uniform_generator(MIN, MAX));

		benchmark::request("th::max (2)",
			[MIN, MAX](real x) { return max(x, MAX); }, uniform_generator(MIN, MAX));

		benchmark::request("th::min (1)",
			[MIN, MAX](real x) { return min(MIN, x); }, uniform_generator(MIN, MAX));

		benchmark::request("th::min (2)",
			[MIN, MAX](real x) { return min(x, MAX); }, uniform_generator(MIN, MAX));

		benchmark::request("th::clamp (1)",
			[MIN, MAX](real x) { return clamp(x, MIN, MAX); }, uniform_generator(MIN, MAX));

		benchmark::request("th::clamp (2)",
			[MIN, MAX](real x) { return clamp(x, 0, 1); }, uniform_generator(MIN, MAX));

		BENCHMARK(th::ln, 0, MAX);
		BENCHMARK(th::log2, 0, MAX);
		BENCHMARK(th::log10, 0, MAX);

		benchmark::request("th::pow (1)",
			[](real x) { return th::pow(1.1, (int) x); }, uniform_generator(-100, 100));

		benchmark::request("th::pow (2)",
			[](real x) { return th::pow(1.1, (int) -x); }, uniform_generator(-100, 100));

		benchmark::request("th::root",
			[](real x) { return th::root(x, 10); }, uniform_generator(MIN, MAX), 100000, 5);

		BENCHMARK(th::fract, 0, 20);
		BENCHMARK(th::exp, -100, 10);

		benchmark::request("th::powf (1)",
			[](real x) { return th::powf(x, 10); }, uniform_generator(MIN, MAX));

		benchmark::request("th::powf (2)",
			[](real x) { return th::powf(x, -10); }, uniform_generator(MIN, MAX));

		BENCHMARK(th::sin, MIN, MAX);
		BENCHMARK(th::cos, MIN, MAX);
		BENCHMARK(th::tan, MIN, MAX);
		BENCHMARK(th::cot, MIN, MAX);

		BENCHMARK(th::atan, MIN, MAX);
		BENCHMARK(th::asin, MIN, MAX);
		BENCHMARK(th::acos, MIN, MAX);

		BENCHMARK(th::sinh, -50, 50);
		BENCHMARK(th::cosh, -50, 50);
		BENCHMARK(th::tanh, -50, 50);
		BENCHMARK(th::coth, -50, 50);
		BENCHMARK(th::sigmoid, -50, 50);
		BENCHMARK(th::sinc, MIN, MAX);
		BENCHMARK(th::heaviside, MIN, MAX);
		BENCHMARK(th::radians, MIN, MAX);
		BENCHMARK(th::degrees, MIN, MAX);

		benchmark::request("th::binomial_coeff",
			[](real x) { return binomial_coeff<uint32_t>(10, x); }, uniform_generator(0, 9));

	benchmark::terminate();
}
