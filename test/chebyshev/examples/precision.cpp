
///
/// @file precision.cpp Example program for precision testing.
///

#include "chebyshev.h"
using namespace ch;


double f(double x) {
	return x * std::sqrt(x);
}


double f_a(double x) {
	return x * std::sqrt(x + 1E-12);
}


int main(int argc, char const *argv[]) {

	// Setup the precision testing environment
	prec::setup("chebyshev", argc, argv);

		// Estimate errors on f_a on [0, 1000]
		prec::estimate("f_a", f_a, f, prec::interval(0, 100));

		// Check that two values are equal
		// up to a tolerance
		prec::equals("f_a", f_a(1), 1, 0.2);

	// Stop precision testing
	prec::terminate();
}
