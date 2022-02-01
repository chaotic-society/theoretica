#ifndef UROBORO_APPROX_H
#define UROBORO_APPROX_H

#include "./vec_buff.h"
#include "./function.h"
#include "./calculus/derivation.h"


namespace uroboro {


	// Root approximation


	// Approximate a root of an arbitrary function using bisection
	// inside a compact interval [a, b]
	inline real approx_root_bisection(real_function f, real a, real b) {

		real x_avg = a;

		real x_min = a;
		real x_max = b;

		int i = 0;

		while(x_max - x_min > BISECTION_APPROX_TOL && i < MAX_BISECTION_ITER) {

			x_avg = (x_max + x_min) / 2.0;

			if(f(x_avg) * f(x_min) > 0)
				x_min = x_avg;
			else
				x_max = x_avg;

			++i;
		}

		return x_avg;
	}


	// Approximate a root of an arbitrary function using Newthon's method
	inline real approx_root_newton(real_function f, real_function Df, real guess = 0) {


		real x = guess;
		int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter < MAX_NEWTON_ITER) {
			x = x - (f(x) / Df(x));
			iter++;
		}

		return x;
	}


	// Approximate a root of a polynomial using Newthon's method
	inline real approx_polyn_root_newton(polynomial p, real guess = 0) {

		real x = guess;
		polynomial Dp = differentiate_polynomial(p);
		int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter < MAX_NEWTON_ITER) {
			x = x - (p(x) / Dp(x));
			iter++;
		}

		return x;
	}


	// Approximate a root of an arbitrary function using Steffensen's method
	inline real approx_root_steffensen(real_function f, real_function Df, real guess = 0) {


		real x = guess;
		int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter < MAX_STEFFENSEN_ITER) {
			x = x - (f(x) / ((f(x + f(x)) / f(x)) - 1));
			iter++;
		}

		return x;
	}


	// Approximate a root of a polynomial using Steffensen's method
	inline real approx_polyn_root_steffensen(polynomial p, real guess = 0) {

		real x = guess;
		int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter < MAX_STEFFENSEN_ITER) {
			x = x - (p(x) / ((p(x + p(x)) / p(x)) - 1));
			iter++;
		}

		return x;
	}


	// Approximate a root of an arbitrary function using Chebyshev's method
	inline real approx_root_chebyshev(real_function f, real_function Df,
		real_function D2f, real guess = 0) {


		real x = guess;
		int iter = 0;

		while(abs(f(x)) > ROOT_APPROX_TOL && iter < MAX_CHEBYSHEV_ITER) {
			x = x - (f(x) / Df(x)) - square((f(x) / Df(x))) * (Df(x) / (D2f(x) * 2));
			iter++;
		}

		return x;
	}


	// Approximate a root of a polynomial using Chebyshev's method
	inline real approx_polyn_root_chebyshev(polynomial p, real guess = 0) {

		real x = guess;
		polynomial Dp = differentiate_polynomial(p);
		polynomial D2p = differentiate_polynomial(p);
		int iter = 0;

		while(abs(p(x)) > ROOT_APPROX_TOL && iter < MAX_CHEBYSHEV_ITER) {
			x = x - (p(x) / Dp(x)) - square((p(x) / Dp(x))) * (Dp(x) / (D2p(x) * 2));
			iter++;
		}

		return x;
	}


	// Max & Min approximation


	// Approximate a function maximum given the function and the first
	// two derivatives using Newton-Raphson
	inline real approx_max_newton(
		real_function f, real_function Df, real_function D2f,
		real guess = 0, real dx = 0.01) {

		real z = approx_root_newton(Df, D2f, guess);

		if(D2f(z) < 0)
			return z;

		// TO-DO throw...
		return 0;
	}


	// Approximate a function minimum given the function and the first
	// two derivatives using Newton-Raphson
	inline real approx_min_newton(
		real_function f, real_function Df,
		real_function D2f, real guess = 0) {

		real z = approx_root_newton(Df, D2f, guess);

		if(D2f(z) > 0)
			return z;

		// TO-DO throw...
		return 0;
	}


	// Approximate a function maximum inside an interval given
	// the function and its first derivative using bisection (on the derivative)
	inline real approx_max_bisection(
		real_function f, real_function Df,
		real a, real b) {

		real z = approx_root_bisection(Df, a, b);

		if(approx_derivative(Df, z) < 0)
			return z;

		// TO-DO throw...
		return 0;
	}


	// Approximate a function minimum inside an interval given
	// the function and its first derivative using bisection (on the derivative)
	inline real approx_min_bisection(real_function f, real_function Df, real a, real b) {

		real z = approx_root_bisection(Df, a, b);

		if(approx_derivative(Df, z) > 0)
			return z;

		// TO-DO throw...
		return 0;
	}


}

#endif
