
#include <cmath>
#include <ctime>
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;



int main(int argc, char const *argv[]) {

	const real MAX = 1E+6;
	PRNG g = PRNG(time(nullptr));
	pdf_sampler gauss = pdf_sampler::gaussian(0, MAX, g);


	prec::state.outputFolder = "test/";
	
	prec::setup("complex");


		// Re() and Im()
		{
			real x = gauss();
			real y = gauss();

			prec::equals("complex::Re()", complex<real>(x, y).Re(), x);
			prec::equals("complex::Re()", complex<real>(x, 0).Re(), x);
			prec::equals("complex::Re()", complex<real>(0).Re(), 0);

			prec::equals("complex::Im()", complex<real>(x, y).Im(), y);
			prec::equals("complex::Im()", complex<real>(0, y).Im(), y);
			prec::equals("complex::Im()", complex<real>(0).Im(), 0);
		}


		// operator+
		{
			real x = gauss();
			real y = gauss();

			prec::equals(
				"complex::operator+",
				(complex<real>(x) + complex<real>(y)).Re(),
				x + y
			);

			prec::equals(
				"complex::operator+",
				(complex<real>(x) + y).Re(),
				x + y
			);

			prec::equals(
				"complex::operator+",
				(x + complex<real>(y)).Re(),
				x + y
			);
		}

	prec::terminate();
}
