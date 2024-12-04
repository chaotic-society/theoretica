
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;
using namespace autodiff;


template<typename Number>
Number f_1D(Number x) {

	return th::exp(x - x * x) + th::sin(x * x) / th::ln(x);
}


template<typename Number>
Number g_1D(Number x) {

	return th::exp(th::sin(x * x)) - 1;
}


dvec2 f_2D(dvec2 v) {

	const dreal2 x = v[0];
	const dreal2 y = v[1];

	return {
		exp(x) - y,
		x * y - exp(x)
	};
}


dvec2 g_2D(dvec2 v) {

    return {
        th::sin(v[0]) - v[1] - 0.5,
        v[0] - v[1]
    };
}


int main(int argc, char const *argv[]) {
	
	prec::setup("optimization");

		prec::settings.outputFiles = { "test/prec/prec_optimization.csv" };

	{
		const real f1_root = 0.690389757422;

		prec::equals(
			"root_bisect (1)",
			root_bisect(f_1D<real>, 0.5, 0.8, 1E-12),
			f1_root
		);
		
		prec::equals("root_newton (1)", root_newton(f_1D, 0.5), f1_root);
		prec::equals("root_halley (1)", root_halley(f_1D<dual2>, 0.5), f1_root);
		prec::equals("root_chebyshev (1)", root_chebyshev(f_1D, 0.5), f1_root);
	}


	{
		const real g1_root = 1.772453850906;

		prec::equals(
			"root_bisect (2)",
			root_bisect(g_1D<real>, 1.5, 2, 1E-12),
			g1_root
		);

		prec::equals("root_newton (2)", root_newton(g_1D, 1.5), g1_root);
		prec::equals("root_halley (2)", root_halley(g_1D<dual2>, 1.5), g1_root);
		prec::equals("root_steffensen (2)", root_steffensen(g_1D<real>, 1), g1_root);
		prec::equals("root_chebyshev (2)", root_chebyshev(g_1D, 1.5), g1_root);
	}


	{
		const auto res = multiroot_newton<2>(f_2D, {1, 1}, 10E-10);
		const auto residual = res - vec2({1, E});
		prec::equals("multiroot_newton (1)", algebra::norm(residual), 0.0);
	}


	{
		const auto res = multiroot_newton<2>(g_2D, {1, 1}, 10E-10);
		const auto residual = res - vec2({-1.49730038909589, -1.49730038909589});
		prec::equals("multiroot_newton (2)", algebra::norm(residual), 0.0);
	}


	prec::terminate();
}
