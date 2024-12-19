
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;
using namespace autodiff;


template<typename Number>
Number f1(Number x) {

	return th::exp(x - x * x) + th::sin(x * x) / th::ln(x);
}


real df1(real x) {

	return th::exp(x - x * x) * (1 - 2 * x)
		+ (2 * x * x * th::cos(x * x) * th::ln(x) - th::sin(x * x))
			/ (x * square(th::ln(x)));
}


template<typename Number>
Number g1(Number x) {

	return th::exp(th::sin(x * x)) - 1;
}


real dg1(real x) {
	return 2 * x * th::exp(th::sin(x * x)) * th::cos(x * x);
}

real d2g1(real x) {
	return 2 * th::exp(th::sin(x * x))
		* (-2 * x * x * th::sin(x * x)
			+ 2 * x * x * square(th::cos(x * x)) + th::cos(x * x));
}


dvec2 f2(dvec2 v) {

	const dreal2 x = v[0];
	const dreal2 y = v[1];

	return {
		exp(x) - y,
		x * y - exp(x)
	};
}


dvec2 g2(dvec2 v) {

    return {
        th::sin(v[0]) - v[1] - 0.5,
        v[0] - v[1]
    };
}


int main(int argc, char const *argv[]) {
	
	prec::setup("optimization");

		prec::settings.outputFiles = { "test/prec/prec_optimization.csv" };
		prec::settings.verboseOutput = true;

	{
		// root = 0.690389757422;

		prec::equals(
			"root_bisect (1)",
			f1(root_bisect(f1<real>, 0.6, 0.7, 1E-12)),
			0.0
		);

		prec::equals(
			"root_itp (1)",
			f1(root_itp(f1<real>, 0.6, 0.7, 1E-12)),
			0.0
		);
		
		prec::equals(
			"root_newton (1)",
			f1(root_newton(f1<dual>, 0.5)), 0.0
		);

		prec::equals(
			"root_newton (1)",
			f1(root_newton(f1<real>, df1, 0.5)), 0.0
		);

		prec::equals(
			"root_halley (1)",
			f1(root_halley(f1<dual2>, 0.5)),
			0.0
		);

		prec::equals(
			"root_chebyshev (1)",
			f1(root_chebyshev(f1<dual2>, 0.7)), 0.0
		);

		prec::equals(
			"root_ostrowski (1)",
			f1(root_ostrowski(f1<real>, df1, 0.7)), 0.0
		);

		prec::equals(
			"root_jarrat (1)",
			f1(root_jarrat(f1<real>, df1, 0.7)), 0.0
		);
	}


	{
		// root = 1.772453850906;

		prec::equals(
			"root_bisect (2)",
			g1(root_bisect(g1<real>, 1.5, 2, 1E-12)),
			0.0
		);

		prec::equals(
			"root_itp (2)",
			g1(root_itp(g1<real>, 1.5, 2, 1E-12)),
			0.0
		);

		prec::equals(
			"root_newton (2)",
			g1(root_newton(g1<dual>, 1.5)), 0.0
		);

		prec::equals(
			"root_newton (2)",
			g1(root_newton(g1<real>, dg1, 1.5)), 0.0
		);
		
		prec::equals(
			"root_halley (2)",
			g1(root_halley(g1<dual2>, 1.5)), 0.0
		);

		prec::equals(
			"root_halley (2)",
			g1(root_halley(g1<real>, dg1, d2g1, 1.5)), 0.0
		);
		
		prec::equals(
			"root_chebyshev (2)",
			g1(root_chebyshev(g1, 1.5)), 0.0
		);

		prec::equals(
			"root_chebyshev (2)",
			g1(root_chebyshev(g1, dg1, d2g1, 1.5)), 0.0
		);

		prec::equals(
			"root_ostrowski (2)",
			g1(root_ostrowski(g1, dg1, 1.5)), 0.0
		);

		prec::equals(
			"root_jarrat (2)",
			g1(root_jarrat(g1, dg1, 1.5)), 0.0
		);
	}



	{
		const auto res = multiroot_newton<2>(f2, {1, 1}, 10E-10);
		const auto residual = res - vec2({1, E});
		prec::equals("multiroot_newton (1)", algebra::norm(residual), 0.0);
	}


	{
		const auto res = multiroot_newton<2>(g2, {1, 1}, 10E-10);
		const auto residual = res - vec2({-1.49730038909589, -1.49730038909589});
		prec::equals("multiroot_newton (2)", algebra::norm(residual), 0.0);
	}


	prec::terminate();
}
