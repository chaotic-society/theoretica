
#include "theoretica.h"
#include <cmath>
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;



int main(int argc, char const *argv[]) {

	prec::settings.outputFiles = { "test/prec/prec_interpolation.csv" };
	
	prec::setup("interpolation");

		// polynomial.h

		// splines.h

	prec::terminate();
}
