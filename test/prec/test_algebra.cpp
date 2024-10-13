
/// @file test_algebra.cpp Test cases for linear algebra

#include "theoretica.h"
#include "chebyshev.h"
#include <ctime>

using namespace theoretica;
using namespace chebyshev;


int main(int argc, char const *argv[]) {

	prec::settings.outputFiles = { "test/prec/test_algebra.csv" };

	prec::setup("algebra", argc, argv);

		// algebra.h

		// mat.h

		// vec.h

		// distance.h
	
	prec::terminate();
}
