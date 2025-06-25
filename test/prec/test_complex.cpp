
#include <cmath>
#include <ctime>
#include "theoretica.h"
#include "prec.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	// Variance for random number generation
	const real VARIANCE = 1E+6;

	auto ctx = prec::make_context("complex", argc, argv);
	ctx.output->settings.outputFiles = { "test/prec/prec_complex.csv" };
	random::random_source rnd = ctx.random->get_rnd();

	// Re() and Im()
	{
		real x = rnd.gaussian(0, VARIANCE);
		real y = rnd.gaussian(0, VARIANCE);

		ctx.equals("complex::Re()", complex<real>(x, y).Re(), x);
		ctx.equals("complex::Re()", complex<real>(x, 0).Re(), x);
		ctx.equals("complex::Re()", complex<real>(0).Re(), 0);

		ctx.equals("complex::Im()", complex<real>(x, y).Im(), y);
		ctx.equals("complex::Im()", complex<real>(0, y).Im(), y);
		ctx.equals("complex::Im()", complex<real>(0).Im(), 0);
	}


	// operator+
	{
		real x = rnd.gaussian(0, VARIANCE);
		real y = rnd.gaussian(0, VARIANCE);

		ctx.equals(
			"complex::operator+",
			(complex<real>(x) + complex<real>(y)).Re(),
			x + y
		);

		ctx.equals(
			"complex::operator+",
			(complex<real>(x) + y).Re(),
			x + y
		);

		ctx.equals(
			"complex::operator+",
			(x + complex<real>(y)).Re(),
			x + y
		);
	}
}
