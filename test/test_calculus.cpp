
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


real f(real x) {
	return th::cos(square(x)) / th::exp(-square(x)) / ln(1 / square(x));
}


real Df(real x) {
	return (2 * th::exp(square(x)) * ((square(x) * ln(1 / square(x)) + 1)
				* th::cos(square(x)) - square(x) * ln(1 / square(x)) * th::sin(square(x))))
					/ (x * square(ln(1 / square(x))));
}


real g(real x) {
	return x * ln(1 / square(x));
}


real G(real x) {
	return 0.5 * square(x) * (ln(1 / square(x)) + 1);
}


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("calculus");

		// Compare the numerical derivative to the analytical derivative

		prec::estimate("deriv_ridders2",
			[](real x) {
				return deriv_ridders2(f, x, 0.00001);
			}, Df,
			interval(0.001, 0.5), prec::state.defaultTolerance, false, 100000, prec::fail_on_rel_err);


		prec::estimate("deriv_ridders",
			[](real x) {
				return deriv_ridders(f, x, 0.00001);
			}, Df,
			interval(0.001, 0.5), prec::state.defaultTolerance, false, 100000, prec::fail_on_rel_err);


		prec::estimate("integral_simpson",
			[](real x) {
				return integral_simpson(g, 1, x);
			},
			[](real x) {
				return G(x) - G(1); // g and G are undefined at 0
			},
			interval(0.1, 3), prec::state.defaultTolerance, false, 100000, prec::fail_on_rel_err);


	prec::terminate();
}
