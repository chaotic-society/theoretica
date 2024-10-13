
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	benchmark::setup("algebra", argc, argv);

		benchmark::settings.outputFiles = { "test/benchmark/benchmark_algebra.h" };
		

	benchmark::terminate();
}
