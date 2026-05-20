
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


real noisy(real x) {
    return x + th::sin(100 * x) / 10;
}


template<typename Number>
Number min1d(Number x) {
	return square(x) + 0.1 * square(th::sin(x));
}

real dmin1d(real x) {
	return 2 * x + 0.1 * th::sin(2 * x);
}

real d2min1d(real x) {
	return 2 + 0.2 * th::cos(2 * x);
}


template<typename Number>
Number max1d(Number x) {
	return -square(x) - 0.1 * square(th::sin(x));
}

real dmax1d(real x) {
	return -2 * x - 0.1 * th::sin(2 * x);
}

real d2max1d(real x) {
	return -2 - 0.2 * th::cos(2 * x);
}


template<typename Number>
Number min2d(vec<Number, 2> v) {

	const Number x = v[0] - 1;
	const Number y = v[1] - 2;
	return square(x) + square(y) + 0.1 * square(th::sin(x + y));
}


template<typename Number>
Number max2d(vec<Number, 2> v) {
	return -min2d(v);
}


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("optimization", argc, argv);
	ctx.settings.outputFiles = { "test/prec/prec_optimization.csv" };
	auto vec2_opt = prec::equation_options<vec2>(
		ctx.settings.defaultTolerance,
		prec::distance::euclidean<vec2>
	);


	// roots.h
	
	{
		// root = 0.690389757422;

		ctx.equals(
			"root_bisect (1)",
			root_bisect(f1<real>, 0.6, 0.7).residual,
			0.0
		);

		ctx.equals(
			"root_itp (1)",
			root_itp(f1<real>, 0.6, 0.7).residual,
			0.0
		);

		ctx.equals(
			"root (1)",
			root(f1<real>, 0.6, 0.7).residual,
			0.0
		);
		
		ctx.equals(
			"root_newton (1)",
			root_newton(f1<dual>, 0.5).residual,
			0.0
		);

		ctx.equals(
			"root_newton (1)",
			root_newton(f1<real>, df1, 0.5).residual,
			0.0
		);

		ctx.equals(
			"root_halley (1)",
			root_halley(f1<dual2>, 0.5).residual,
			0.0
		);

		ctx.equals(
			"root_chebyshev (1)",
			root_chebyshev(f1<dual2>, 0.7).residual, 0.0
		);

		ctx.equals(
			"root_ostrowski (1)",
			root_ostrowski(f1<real>, df1, 0.7).residual, 0.0
		);

		ctx.equals(
			"root_jarrat (1)",
			root_jarrat(f1<real>, df1, 0.7).residual, 0.0
		);
	}


	{
		// root = 1.772453850906;

		ctx.equals(
			"root_newton (2)",
			root_newton(g1<dual>, 1.5).residual, 0.0
		);

		ctx.equals(
			"root_newton (2)",
			root_newton(g1<real>, dg1, 1.5).residual, 0.0
		);
		
		ctx.equals(
			"root_halley (2)",
			root_halley(g1<dual2>, 1.5).residual, 0.0
		);

		ctx.equals(
			"root_halley (2)",
			root_halley(g1<real>, dg1, d2g1, 1.5).residual, 0.0
		);
		
		ctx.equals(
			"root_chebyshev (2)",
			root_chebyshev(g1, 1.5).residual, 0.0
		);

		ctx.equals(
			"root_chebyshev (2)",
			root_chebyshev(g1, dg1, d2g1, 1.5).residual, 0.0
		);

		ctx.equals(
			"root_ostrowski (2)",
			root_ostrowski(g1, dg1, 1.5).residual, 0.0
		);

		ctx.equals(
			"root_jarrat (2)",
			root_jarrat(g1, dg1, 1.5).residual, 0.0
		);
	}


	{
		// root = 0.0

		ctx.equals(
			"root_bisect (2)",
			root_bisect(noisy, -0.031, +0.034).residual,
			0.0
		);

		ctx.equals(
			"root_itp (2)",
			root_itp(noisy, -0.032, +0.034).residual,
			0.0
		);

		ctx.equals(
			"root (2)",
			root(noisy, -0.032, +0.033).residual,
			0.0
		);
	}


	// multi_roots.h
	
	{
		ctx.equals(
			"multiroot_newton (1)",
			multiroot_newton(f2, vec2({1, 1}), 10E-10).residual,
			0.0
		);
	}


	{
		ctx.equals(
			"multiroot_newton (2)",
			multiroot_newton(g2, vec2({1, 1}), 10E-10).residual,
			0.0
		);
	}


	// extrema.h

	{
		auto result = maximize_golden(max1d<real>, 0.0, 4.0);

		ctx.equals(
			"maximize_golden (1)",
			result.value,
			0.0,
			1E-07
		);

		ctx.equals(
			"maximize_golden (1) converged",
			result.converged(),
			true
		);
	}


	{
		auto result = minimize_golden(min1d<real>, 0.0, 4.0);

		ctx.equals(
			"minimize_golden (1)",
			result.value,
			0.0,
			1E-07
		);

		ctx.equals(
			"minimize_golden (1) converged",
			result.converged(),
			true
		);
	}


	{
		auto result = maximize_newton(max1d<real>, dmax1d, d2max1d, 0.0);

		ctx.equals(
			"maximize_newton (1)",
			result.value,
			0.0
		);

		ctx.equals(
			"maximize_newton (1) converged",
			result.converged(),
			true
		);
	}


	{
		auto result = minimize_newton(min1d<real>, dmin1d, d2min1d, 0.0);

		ctx.equals(
			"minimize_newton (1)",
			result.value,
			0.0
		);

		ctx.equals(
			"minimize_newton (1) converged",
			result.converged(),
			true
		);
	}


	{
		auto result = maximize_bisection(max1d<real>, dmax1d, -2.0, 2.0);

		ctx.equals(
			"maximize_bisection (1)",
			result.value,
			0.0,
			1E-07
		);

		ctx.equals(
			"maximize_bisection (1) converged",
			result.converged(),
			true
		);
	}


	{
		auto result = minimize_bisect(min1d<real>, dmin1d, -0.5, 2.0);

		ctx.equals(
			"minimize_bisect (1)",
			result.value,
			0.0,
			1E-07
		);

		ctx.equals(
			"minimize_bisect (1) converged",
			result.converged(),
			true
		);
	}


	// multi_extrema.h

	{
		auto result = multi_minimize_grad(min2d<dreal2>, vec2({0, 0}), 0.1);

		ctx.equals(
			"multi_minimize_grad (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}


	{
		auto result = multi_maximize_grad(max2d<dreal2>, vec2({0, 0}), 0.1);

		ctx.equals(
			"multi_maximize_grad (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}


	{
		auto result = multi_minimize_lingrad(min2d<dreal2>, vec2({0, 0}));

		ctx.equals(
			"multi_minimize_lingrad (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}


	{
		auto result = multi_maximize_lingrad(max2d<dreal2>, vec2({0, 0}));

		ctx.equals(
			"multi_maximize_lingrad (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}


	{
		auto result = multi_minimize(min2d<dreal2>, vec2({0, 0}));

		ctx.equals(
			"multi_minimize (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}


	{
		auto result = multi_maximize(max2d<dreal2>, vec2({0, 0}));

		ctx.equals(
			"multi_maximize (1)",
			result,
			vec2({1, 2}),
			vec2_opt
		);
	}
}
