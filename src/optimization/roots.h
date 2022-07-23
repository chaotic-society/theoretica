
///
/// @file roots.h Root approximation of real functions
///

#ifndef THEORETICA_ROOTS_H
#define THEORETICA_ROOTS_H

#include "../core/function.h"
#include "../calculus/derivation.h"
#include "../autodiff/dual.h"
#include "../autodiff/dual2.h"


namespace theoretica {


	/// Approximate a root of an arbitrary function using bisection
	/// inside a compact interval [a, b] where f(a) * f(b) < 0
	inline real approx_root_bisection(real_function f, real a, real b) {

		if(f(a) * f(b) >= 0) {
			TH_MATH_ERROR("approx_root_bisection", f(a) * f(b), INVALID_ARGUMENT);
			return nan();
		}

		real x_avg = a;

		real x_min = a;
		real x_max = b;

		unsigned int iter = 0;

		while(x_max - x_min > BISECTION_APPROX_TOL && iter <= MAX_BISECTION_ITER) {

			x_avg = (x_max + x_min) / 2.0;

			if(f(x_avg) * f(x_min) > 0)
				x_min = x_avg;
			else
				x_max = x_avg;

			++iter;
		}

		if(iter > MAX_BISECTION_ITER) {
			TH_MATH_ERROR("approx_root_bisection", x_avg, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x_avg;
	}


	/// Approximate a root of an arbitrary function using Newthon's method
	inline real approx_root_newton(real_function f, real_function Df, real guess = 0) {


		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > NEWTON_RAPHSON_TOL && iter <= MAX_NEWTON_ITER) {
			x = x - (f(x) / Df(x));
			iter++;
		}

		if(iter > MAX_NEWTON_ITER) {
			TH_MATH_ERROR("approx_root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Newthon's method,
	/// computing the derivative using automatic differentiation
	inline real approx_root_newton(dual(*f)(dual), real guess = 0) {


		real x = guess;
		unsigned int iter = 0;

		dual s = f(dual(x, 1));

		while(abs(s.Re()) > NEWTON_RAPHSON_TOL && iter <= MAX_NEWTON_ITER) {

			s = f(dual(x, 1));

			x = x - (s.Re() / s.Dual());
			iter++;
		}

		if(iter > MAX_NEWTON_ITER) {
			TH_MATH_ERROR("approx_root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Newthon's method
	inline real approx_polyn_root_newton(polynomial<real> p, real guess = 0) {

		real x = guess;
		polynomial<> Dp = deriv_polynomial(p);
		unsigned int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter <= MAX_NEWTON_ITER) {
			x = x - (p(x) / Dp(x));
			iter++;
		}

		if(iter > MAX_NEWTON_ITER) {
			TH_MATH_ERROR("approx_polyn_root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Halley's method
	inline real approx_root_halley(real_function f, real_function Df,
		real_function D2f, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter <= MAX_HALLEY_ITER) {
			x = x - (2 * f(x) * Df(x)) / (2 * Df(x) - f(x) * D2f(x));
			iter++;
		}

		if(iter > MAX_HALLEY_ITER) {
			TH_MATH_ERROR("approx_root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Halley's method
	inline real approx_root_halley(dual2(*f)(dual2), real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		dual2 s = f(dual2(x, 1, 0));

		while(abs(s.Re()) > ROOT_APPROX_TOL && iter <= MAX_HALLEY_ITER) {

			s = f(dual2(x, 1, 0));

			x = x - (2 * s.Re() * s.Dual1())
						/ (2 * s.Dual1() - s.Re() * s.Dual2());
			iter++;
		}

		if(iter > MAX_HALLEY_ITER) {
			TH_MATH_ERROR("approx_root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Halley's method
	inline real approx_polyn_root_halley(polynomial<real> p, real guess = 0) {

		polynomial<> Dp = deriv_polynomial(p);
		polynomial<> D2p = deriv_polynomial(Dp);
		real x = guess;
		unsigned int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter <= MAX_HALLEY_ITER) {
			x = x - (2 * p(x) * Dp(x)) / (2 * Dp(x) - p(x) * D2p(x));
			iter++;
		}

		if(iter > MAX_HALLEY_ITER) {
			TH_MATH_ERROR("approx_polyn_root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Steffensen's method
	inline real approx_root_steffensen(real_function f, real guess = 0) {


		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter <= MAX_STEFFENSEN_ITER) {
			x = x - (f(x) / ((f(x + f(x)) / f(x)) - 1));
			iter++;
		}

		if(iter > MAX_STEFFENSEN_ITER) {
			TH_MATH_ERROR("approx_root_steffensen", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Steffensen's method
	inline real approx_polyn_root_steffensen(polynomial<real> p, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter <= MAX_STEFFENSEN_ITER) {
			x = x - (p(x) / ((p(x + p(x)) / p(x)) - 1));
			iter++;
		}

		if(iter > MAX_STEFFENSEN_ITER) {
			TH_MATH_ERROR("approx_polyn_root_steffensen", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Chebyshev's method
	inline real approx_root_chebyshev(real_function f, real_function Df,
		real_function D2f, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter <= MAX_CHEBYSHEV_ITER) {
			x = x - (f(x) / Df(x)) - square((f(x) / Df(x))) * (Df(x) / (D2f(x) * 2));
			iter++;
		}

		if(iter > MAX_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("approx_root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Chebyshev's method
	inline real approx_root_chebyshev(dual2(*f)(dual2), real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		dual2 s = f(dual2(x, 1, 0));

		while(abs(s.Re()) > ROOT_APPROX_TOL && iter <= MAX_CHEBYSHEV_ITER) {

			s = f(dual2(x, 1, 0));

			x = x - (s.Re() / s.Dual1())
						- square((s.Re() / s.Dual1()))
							* (s.Dual1() / (s.Dual2() * 2));
			iter++;
		}

		if(iter > MAX_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("approx_root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Chebyshev's method
	inline real approx_polyn_root_chebyshev(polynomial<real> p, real guess = 0) {

		real x = guess;
		polynomial<> Dp = deriv_polynomial(p);
		polynomial<> D2p = deriv_polynomial(p);
		unsigned int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter <= MAX_CHEBYSHEV_ITER) {
			x = x - (p(x) / Dp(x)) - square((p(x) / Dp(x))) * (Dp(x) / (D2p(x) * 2));
			iter++;
		}

		if(iter > MAX_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("approx_polyn_root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}

}

#endif
