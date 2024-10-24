
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


// Multivariate root search
// {1, e} is a root
dvec2 f_2D(dvec2 v) {

	const dreal2 x = v[0];
	const dreal2 y = v[1];

	return {
		exp(x) - y,
		x * y - exp(x)
	};
}


int main(int argc, char const *argv[]) {
	
	prec::setup("optimization");

		prec::settings.outputFiles = { "test/prec/prec_optimization.csv" };

		// const real f1_root = 0.690389757422;
		// const real g1_root = 1.772453850906;

		// prec::equals(
		// 	"root_bisection",
		// 	root_bisection(f_1D<real>, 0.5, 0.8, 1E-12),
		// 	f1_root
		// );
		
		// prec::equals("root_newton", root_newton(f_1D, 0.5), f1_root);
		// prec::equals("root_chebyshev", root_chebyshev(f_1D, 0.5), f1_root);


		// prec::equals(
		// 	"root_bisection",
		// 	root_bisection(g_1D<real>, 1.5, 2, 1E-12),
		// 	g1_root
		// );

		// prec::equals("root_newton", root_newton(g_1D, 1.5), g1_root);
		// prec::equals("root_steffensen", root_steffensen(g_1D<real>, 1), g1_root);
		// prec::equals("root_chebyshev", root_chebyshev(g_1D, 1.5), g1_root);


		// const auto r = multiroot_newton<2>(f_2D, {2, 2}, 10E-10);
		// prec::equals("multiroot_newton (1)", r[0], 1);
		// prec::equals("multiroot_newton (2)", r[1], E);


	prec::terminate();
}
