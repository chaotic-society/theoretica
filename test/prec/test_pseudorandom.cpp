
#include "theoretica.h"
#include <cmath>
#include "chebyshev.h"
#include <ctime>

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("pseudorandom", argc, argv);
	ctx.settings.outputFiles = { "test/prec/prec_pseudorandom.csv" };
	
	// pseudorandom.h

	// prng.h

	// montecarlo.h

	// sampling.h
}
