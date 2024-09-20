
#include "theoretica.h"
#include <cmath>
#include "chebyshev.h"

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


vec2 diff_eq(real t, vec2 v) {

	return {
		v[1],
		-v[0]
	};
}


vec2 sho(real t) {

	return {
		std::sin(t),
		std::cos(t)
	};
}


int main(int argc, char const *argv[]) {
	
	prec::setup("calculus");

		output::state.outputFiles = { "test/prec_calculus.csv" };
		prec::state.estimateColumns = {
			"funcName", "meanErr", "rmsErr", "maxErr", "tolerance", "failed"
		};

		// Compare the numerical derivative to the analytical derivative

		auto deriv_opt = prec::estimate_options<real, real>(
			prec::interval(0.001, 0.5),
			prec::estimator::quadrature1D(),
			10E-04
		);


		prec::estimate("deriv_forward",
			[](real x) { return deriv_forward(f, x, 10E-8); },
			Df, deriv_opt);


		prec::estimate("deriv_backward",
			[](real x) { return deriv_backward(f, x, 10E-8); },
			Df, deriv_opt);


		prec::estimate("deriv_central",
			[](real x) { return deriv_central(f, x, 10E-8); },
			Df, deriv_opt);


		prec::estimate("deriv_ridders2",
			[](real x) { return deriv_ridders2(f, x, 10E-6); },
			Df, deriv_opt);


		prec::estimate("deriv_ridders",
			[](real x) { return deriv_ridders(f, x, 10E-6, 3); },
			Df, deriv_opt);


		// Compare integral quadrature to primitives

		auto integ_opt = prec::estimate_options<real, real>(
			prec::interval(0.1, 3.0),
			prec::estimator::quadrature1D()
		);


		prec::estimate("integral_trapezoid",
			[](real x) { return integral_trapezoid(g, 1, x); },
			[](real x) { return G(x) - G(1); },
			{ prec::interval(0.1, 3.0) },
			1E-04, 1'000, prec::fail::fail_on_max_err(),
			prec::estimator::quadrature1D()
		);


		prec::estimate<real, real>(
			"integral_simpson",
			[](real x) { return integral_simpson(g, 1, x); },
			[](real x) { return G(x) - G(1); },
			integ_opt
		);


		prec::estimate("integral_romberg",
			[](real x) { return integral_romberg(g, 1, x); },
			[](real x) { return G(x) - G(1); },
			integ_opt
		);


		prec::estimate("integral_legendre",
			[](real x) { return integral_legendre(g, 1, x, 16); },
			[](real x) { return G(x) - G(1); },
			integ_opt
		);

	prec::terminate();
}
