
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	auto ctx = benchmark::make_context("algebra", argc, argv);
	ctx.settings.outputFiles = { "test/benchmark/benchmark_algebra.csv" };
	auto rng = ctx.random->get_rnd();
	ctx.settings.multithreading = false;
	
	using Matrix = mat<real>;
	std::vector<std::pair<Matrix, Matrix>> matrices (20);

	ctx.output->info("\tGenerating random matrices...");
	for (size_t i = 0; i < matrices.size(); i++) {

		matrices[i].first.resize(256, 256);
		matrices[i].second.resize(256, 256);
		
		for (auto& x : matrices[i].first)
			x = rng.gaussian(0, 1);

		for (auto& x : matrices[i].second)
			x = rng.gaussian(0, 1);
	}

	ctx.output->info("\tLaunching benchmark...");
	ctx.benchmark(
		"algebra::mat_mul",
		[](const auto& input) {
			algebra::mat_mul(input.first, input.second);
			return input.first(0, 0) + input.second(0, 0);
		},
		matrices, 50
	);
}
