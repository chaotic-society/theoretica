
#include "theoretica.h"
#include <cmath>
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	auto ctx = prec::make_context("interpolation");
	ctx.settings.outputFiles = { "test/prec/prec_interpolation.csv" };

	// polynomial.h

	// splines.h

}
