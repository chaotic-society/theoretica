
///
/// @file test_polynomial.cpp Polynomial class test cases
///

#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	prec::setup("polynomial");

		prec::settings.outputFiles = { "test/prec/prec_polynomial.csv" };

		// polynomial.h

		// orthogonal.h

	prec::terminate();
}
