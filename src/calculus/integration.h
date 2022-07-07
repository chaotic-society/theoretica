
///
/// @file integration.h Integral approximation
///

#ifndef THEORETICA_INTEGRATION_H
#define THEORETICA_INTEGRATION_H

#include "../core/constants.h"
#include "../core/function.h"
#include "../polynomial/polynomial.h"


namespace theoretica {


	// Integrate a polynomial
	polynomial<real> integrate_polynomial(const polynomial<real>& p) {

		polynomial<real> Dp;
		Dp.coeff.push_back(0);

		for (unsigned int i = 0; i < p.size(); ++i)
			Dp.coeff.push_back(p.get(i) / (real) (i + 1));

		return Dp;
	}


	// Approximate the definite integral of an arbitrary function
	// using the midpoint method
	real approx_integral_midpoint(real_function f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {
		
		real dx = (b - a) / steps;
		real res = 0;

		for (unsigned int i = 0; i < steps; ++i)
			res += f(a + (i + 0.5) * dx);

		return res * dx;
	}


	// Approximate the definite integral of an arbitrary function
	// using the trapezoid method
	real approx_integral_trapezoid(real_function f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {
		
		real dx = (b - a) / steps;
		real res = 0;

		res += 0.5 * f(a);

		for (unsigned int i = 1; i < steps; ++i)
			res += f(a + i * dx);

		res += 0.5 * f(b);

		return res * dx;
	}


	// Approximate the definite integral of an arbitrary function
	// using Simpson's method
	real approx_integral_simpson(real_function f, real a, real b,
		unsigned int steps = INTEGRATION_STEPS) {
		
		real dx = (b - a) / (real) steps;
		real res = 0;

		res += f(a) + f(b);

		for (unsigned int i = 1; i < steps; ++i) {

			if(i % 2 == 0)
				res += 2.0 * f(a + i * dx);
			else
				res += 4.0 * f(a + i * dx);
			
		}

		return res * dx / 3.0;
	}


	// Approximate the definite integral of an arbitrary function
	// using Romberg's method accurate to the given order
	real approx_integral_romberg(real_function f, real a, real b, unsigned int order) {

		if(order % 2 != 0) {
			UMATH_ERROR("approx_integral_romberg", order, IMPOSSIBLE_OPERATION);
			return nan();
		}

		unsigned int iter = order / 2;

		real T[iter][iter];

		T[0][0] = (f(a) + f(b)) * (b - a) / 2.0;

		for (unsigned int j = 1; j < iter; ++j) {
			
			// Composite Trapezoidal Rule
			T[j][0] = approx_integral_trapezoid(f, a, b, 1 << j);

			// Richardson extrapolation
			for (unsigned int k = 1; k <= j; ++k) {
				T[j][k] =
					T[j][k - 1]
					+ (T[j][k - 1] - T[j - 1][k - 1]) / ((1 << (2 * k)) - 1);
			}
		}

		// Return the best approximation
		return T[iter - 1][iter - 1];
	}

}


#endif
