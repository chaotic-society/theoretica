
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

		output::state.outputFolder = "test/";
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


		// Integrate the simple harmonic oscillator

		const real h_euler = 1E-4;
		const real h_rk2 = 1E-4;
		const real h_rk4 = 1E-4;

		real time = 1;

		auto ode_opt = prec::estimate_options<real, real>(
			prec::interval(0, time),
			prec::estimator::quadrature1D(),
			1E-03
		);

		ode_state<2> s_euler(0, {0, 1});
		ode_state<2> s_rk2(0, {0, 1});
		ode_state<2> s_rk4(0, {0, 1});

		std::vector<vec2> traj_euler;
		std::vector<vec2> traj_rk2;
		std::vector<vec2> traj_rk4;

		std::vector<real> x(1, 0);
		std::vector<real> y(1, 0);

		traj_euler.push_back({0, s_euler.y[0]});
		traj_rk2.push_back({0, s_rk2.y[0]});
		traj_rk4.push_back({0, s_rk4.y[0]});


		for (int i = 0; i < int(time / h_euler); ++i) {
			s_euler = ode_euler(diff_eq, s_euler, h_euler);
			traj_euler.push_back({s_euler.t, s_euler.y[0]});
			x.push_back(s_euler.t);
			y.push_back(s_euler.y[0]);
		}

		for (int i = 0; i < int(time / h_rk2); ++i) {
			s_rk2 = ode_rk2(diff_eq, s_rk2, h_rk2);
			traj_rk2.push_back({h_rk2 * i, s_rk2.y[0]});
		}

		for (int i = 0; i < int(time / h_rk4); ++i) {
			s_rk4 = ode_rk4(diff_eq, s_rk4, h_rk4);
			traj_rk4.push_back({h_rk4 * i, s_rk4.y[0]});
		}


		spline interp_euler = spline(x, y);
		// spline interp_rk2 = spline(traj_rk2);
		// spline interp_rk4 = spline(traj_rk4);

		prec::estimate(
			"ode_euler",
			[interp_euler](real t) { return interp_euler(t); },
			[](real t) { return sho(t)[0]; },
			ode_opt
		);


		// prec::estimate(
		// 	"ode_rk2",
		// 	[interp_rk2](real t) { return interp_rk2(t); },
		// 	[](real t) { return sho(t)[0]; },
		// 	interval(0, time), 1E-4
		// );


		// prec::estimate(
		// 	"ode_rk4",
		// 	[interp_rk4](real t) { return interp_rk4(t); },
		// 	[](real t) { return sho(t)[0]; },
		// 	interval(0, time), 1E-4
		// );

		// Check for energy conservation

		prec::equals("ode_euler (E)", s_euler.y.sqr_norm(), 1.0, 1E-3);

		prec::equals("ode_rk2 (E)", s_rk2.y.sqr_norm(), 1.0, 1E-4);

		prec::equals("ode_rk4 (E)", s_rk4.y.sqr_norm(), 1.0, 1E-4);


	prec::terminate();
}
