
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	auto ctx = benchmark::make_context("algebra", argc, argv);
	ctx.settings.outputFiles = { "test/benchmark/benchmark_algebra.h" };
	
}
