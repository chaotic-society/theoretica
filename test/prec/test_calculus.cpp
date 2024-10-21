
#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

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


// e^{x} sin(x) e^{-x^2}
real GaussI(real x) {
	return th::exp(x) * th::sin(x);
}

// sqrt(x) e^{-x}
real ExpI(real x) {
	return th::cos(x);
}


template<typename Field>
Field h(Field x) {
	return x * sin(x) - cos(x);
}


vec2 diff_eq(real t, vec2 v) {

	return {
		+v[1],
		-v[0]
	};
}


vec2 sho(real t) {
	return {
		std::sin(t),
		std::cos(t)
	};
}


// Construct a precision estimator for an ODE solution
template<typename Vector = vec<real>>
auto ode_estimator(const ode::ode_solution_t<Vector>& sol) {

	return [sol](
		std::function<Vector(real)> dummy,
		std::function<Vector(real)> exact,
		prec::estimate_options<Vector, real> opt) -> prec::estimate_result {

		real absErr = 0.0;
		real sqrAbsErr = 0.0;
		real maxErr = -inf();

		for (unsigned int i = 0; i < sol.t.size(); ++i) {
			
			const Vector delta = exact(sol.t[i]) - sol.x[i];
			real sum = 0.0;

			for (unsigned int j = 0; j < delta.size(); ++j)
				sum += delta[j] * delta[j];

			const real norm = std::sqrt(sum);
			absErr += norm;
			sqrAbsErr += sum;
			maxErr = std::max(maxErr, norm);
		}

		prec::estimate_result res;
		res.maxErr = maxErr;
		res.absErr = absErr;
		res.meanErr = absErr / sol.t.size();
		res.rmsErr = std::sqrt(sqrAbsErr / sol.t.size());
		return res;
	};
}


long double distance_polyn(const polynomial<real>& p1, const polynomial<real>& p2) {

	const polynomial<real> d = p1 - p2;
	real r = -inf();

	for (size_t i = 0; i < d.size(); ++i)
		r = max(r, std::abs(d[i]));

	return r;
}



int main(int argc, char const *argv[]) {
	
	prec::setup("calculus");

		output::settings.outputFiles = { "test/prec/prec_calculus.csv" };
		output::settings.fieldOptions["name"].columnWidth = 24;
		prec::settings.estimateColumns = {
			"name", "meanErr", "rmsErr", "maxErr", "tolerance", "failed"
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


		prec::estimate("integral_midpoint",
			[](real x) { return integral_midpoint(g, 1, x); },
			[](real x) { return G(x) - G(1); },
			{ prec::interval(0.1, 3.0) },
			1E-04, 1'000, prec::fail::fail_on_max_err(),
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


		prec::estimate("integral_romberg_tol",
			[](real x) { return integral_romberg_tol(g, 1, x); },
			[](real x) { return G(x) - G(1); },
			integ_opt
		);


		prec::estimate("integral_legendre",
			[](real x) { return integral_legendre(g, 1, x, 16); },
			[](real x) { return G(x) - G(1); },
			integ_opt
		);


		prec::equals(
			"integral_hermite",
			integral_hermite(GaussI),
			0.8497596421214707431181
		);


		prec::equals(
			"integral_laguerre",
			integral_laguerre(ExpI),
			0.5
		);


		// Integrate the simple harmonic oscillator

		real tf = 1.0;
		vec2 x0 = {0.0, 1.0};

		auto emptyf = [](real t) -> vec2 { return vec2(); };

		// Custom ODE estimator
		auto opt = prec::estimate_options<vec2, real>();


		// Lower order methods
		opt.tolerance = 1E-04;

		opt.estimator = ode_estimator(
			ode::solve_euler(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_euler", emptyf, sho, opt);


		// Higher order methods
		opt.tolerance = 1E-08;

		opt.estimator = ode_estimator(
			ode::solve_midpoint(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_midpoint", emptyf, sho, opt);


		opt.estimator = ode_estimator(
			ode::solve_heun(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_heun", emptyf, sho, opt);

		opt.estimator = ode_estimator(
			ode::solve_rk2(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_rk2", emptyf, sho, opt);


		opt.estimator = ode_estimator(
			ode::solve_rk4(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_rk4", emptyf, sho, opt);


		opt.estimator = ode_estimator(
			ode::solve_k38(diff_eq, x0, 0.0, tf)
		);
		prec::estimate("ode::solve_k38", emptyf, sho, opt);


		// taylor.h


		auto taylor_opt = prec::equation_options<polynomial<real>>(
			1E-08, distance_polyn
		);

	{
		polynomial<real> evaluated = taylor::expand_linear(h<dual>);
		polynomial<real> expected = {-1, 0};

		prec::equals(
			"taylor::expand_linear",
			evaluated, expected, taylor_opt
		);
	}

	{
		polynomial<real> evaluated = taylor::expand_quadratic(h<dual2>);
		polynomial<real> expected = {-1, 0, 1.5};

		prec::equals(
			"taylor::expand_quadratic",
			evaluated, expected, taylor_opt
		);
	}

	prec::terminate();
}
