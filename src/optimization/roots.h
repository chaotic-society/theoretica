
///
/// @file roots.h Root approximation of real functions
///

#ifndef THEORETICA_ROOTS_H
#define THEORETICA_ROOTS_H

#include "../core/function.h"
#include "../calculus/derivation.h"
#include "../autodiff/dual.h"
#include "../autodiff/dual2.h"
#include "../algebra/vec.h"
#include "../complex/complex.h"


namespace theoretica {


	/// Find candidate intervals for root finding
	///
	/// @param f A function of real variable
	/// @param a The lower extreme of the region of interest
	/// @param b The upper extreme of the region of interest
	/// @param steps The number of subintervals to check (optional)
	template<typename RealFunction>
	inline std::vector<vec2> find_root_intervals(
		RealFunction f, real a, real b, unsigned int steps = 10) {

		std::vector<vec2> res;
		const real dx = (b - a) / (real) steps;

		for (unsigned int i = 0; i < steps; ++i) {
			
			const real x1 = a + i * dx;
			const real x2 = a + (i + 1) * dx;

			if(f(x1) * f(x2) <= 0)
				res.push_back({x1, x2});
		}

		return res;
	}


	/// Approximate a root of an arbitrary function using bisection
	/// inside a compact interval [a, b] where f(a) * f(b) < 0
	///
	/// @param f The real function to search the root of.
	/// @param a The lower extreme of the interval.
	/// @param b The upper extreme of the interval.
	/// @param tolerance The minimum size of the bound to stop the algorithm.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_bisection(
		RealFunction f, real a, real b, real tolerance = OPTIMIZATION_TOL) {

		if(f(a) * f(b) >= 0) {
			TH_MATH_ERROR("root_bisection", f(a) * f(b), INVALID_ARGUMENT);
			return nan();
		}

		real x_avg = a;

		real x_min = a;
		real x_max = b;

		unsigned int iter = 0;

		while((x_max - x_min) > tolerance && iter <= OPTIMIZATION_BISECTION_ITER) {

			x_avg = (x_max + x_min) / 2.0;

			if(f(x_avg) * f(x_min) > 0)
				x_min = x_avg;
			else
				x_max = x_avg;

			++iter;
		}

		if(iter > OPTIMIZATION_BISECTION_ITER) {
			TH_MATH_ERROR("root_bisection", x_avg, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x_avg;
	}


	/// Approximate a root of an arbitrary function using Newton's method
	///
	/// @param f The real function to search the root of.
	/// @param Df The derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_newton(RealFunction f, RealFunction Df, real guess = 0) {


		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_NEWTON_ITER) {
			x = x - (f(x) / Df(x));
			iter++;
		}

		if(iter > OPTIMIZATION_NEWTON_ITER) {
			TH_MATH_ERROR("root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Newton's method,
	/// computing the derivative using automatic differentiation.
	///
	/// @param f The real function to search the root of,
	/// with dual argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline real root_newton(dual(*f)(dual), real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		dual s = f(dual(x, 1));

		while(abs(s.Re()) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_NEWTON_ITER) {

			s = f(dual(x, 1));

			x = x - (s.Re() / s.Dual());
			iter++;
		}

		if(iter > OPTIMIZATION_NEWTON_ITER) {
			TH_MATH_ERROR("root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Newton's method.
	///
	/// @param p The polynomial to search the root of.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the polynomial,
	/// or NaN if the algorithm did not converge.
	inline real root_newton(polynomial<real> p, real guess = 0) {

		real x = guess;
		polynomial<> Dp = deriv(p);
		unsigned int iter = 0;

		while(abs(p(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_NEWTON_ITER) {
			x = x - (p(x) / Dp(x));
			iter++;
		}

		if(iter > OPTIMIZATION_NEWTON_ITER) {
			TH_MATH_ERROR("root_newton", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary complex function
	/// using Newton's method,
	///
	/// @param f The complex function to search the root of
	/// @param df The derivative of the function
	/// @param guess The initial guess (defaults to 0).
	/// @param tolerance The minimum tolerance to stop the algorithm
	/// when reached.
	/// @param max_iter The maximum number of iterations before stopping
	/// the algorithm (defaults to OPTIMIZATION_TOL).
	/// @return The coordinate of the root of the function,
	/// or a complex NaN if the algorithm did not converge.
	inline complex<> root_newton(
		complex<>(*f)(complex<>),
		complex<>(*df)(complex<>),
		complex<> guess = complex<>(0, 0),
		real tolerance = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_TOL) {

		complex<> z = guess;
		unsigned int iter = 0;

		while(abs(f(z)) > tolerance && iter <= max_iter) {
			z = z - (f(z) / df(z));
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_newton", z.Re(), NO_ALGO_CONVERGENCE);
			return nan();
		}

		return z;
	}


	/// Approximate a root of an arbitrary function using Halley's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_halley(RealFunction f, RealFunction Df,
		RealFunction D2f, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_HALLEY_ITER) {
			x = x - (2 * f(x) * Df(x)) / (2 * Df(x) - f(x) * D2f(x));
			iter++;
		}

		if(iter > OPTIMIZATION_HALLEY_ITER) {
			TH_MATH_ERROR("root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Halley's method,
	/// leveraging automatic differentiation to compute the first and second
	/// derivatives of the function.
	///
	/// @param f The real function to search the root of,
	/// with dual2 argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline real root_halley(dual2(*f)(dual2), real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		dual2 s = f(dual2(x, 1.0, 0.0));

		while(abs(s.Re()) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_HALLEY_ITER) {

			s = f(dual2(x, 1, 0));

			x = x - (2 * s.Re() * s.Dual1())
						/ (2 * s.Dual1() - s.Re() * s.Dual2());
			iter++;
		}

		if(iter > OPTIMIZATION_HALLEY_ITER) {
			TH_MATH_ERROR("root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Halley's method.
	///
	/// @param p The polynomial to search the root of.
	/// @guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the polynomial,
	/// or NaN if the algorithm did not converge.
	inline real root_halley(polynomial<real> p, real guess = 0) {

		polynomial<> Dp = deriv(p);
		polynomial<> D2p = deriv(Dp);

		real x = guess;
		unsigned int iter = 0;

		while(abs(p(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_HALLEY_ITER) {
			x = x - (2 * p(x) * Dp(x)) / (2 * Dp(x) - p(x) * D2p(x));
			iter++;
		}

		if(iter > OPTIMIZATION_HALLEY_ITER) {
			TH_MATH_ERROR("root_halley", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Steffensen's method.
	///
	/// @param f The real function to search the root of.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_steffensen(RealFunction f, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_STEFFENSEN_ITER) {

			const real f_x = f(x);
			x = x - (f_x / ((f(x + f_x) / f_x) - 1));
			iter++;
		}

		if(iter > OPTIMIZATION_STEFFENSEN_ITER) {
			TH_MATH_ERROR("root_steffensen", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Steffensen's method.
	///
	/// @param p The polynomial to search the root of.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline real root_steffensen(polynomial<real> p, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(p(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_STEFFENSEN_ITER) {
			x = x - (p(x) / ((p(x + p(x)) / p(x)) - 1));
			iter++;
		}

		if(iter > OPTIMIZATION_STEFFENSEN_ITER) {
			TH_MATH_ERROR("root_steffensen", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Chebyshev's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_chebyshev(RealFunction f, RealFunction Df,
		RealFunction D2f, real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		while(abs(f(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_CHEBYSHEV_ITER) {
			x = x - (f(x) / Df(x)) - square((f(x) / Df(x))) * (Df(x) / (D2f(x) * 2));
			iter++;
		}

		if(iter > OPTIMIZATION_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of an arbitrary function using Chebyshev's method,
	/// by computing the first and second derivatives using automatic differentiation.
	///
	/// @param f The real function to search the root of,
	/// with dual2 argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline real root_chebyshev(dual2(*f)(dual2), real guess = 0) {

		real x = guess;
		unsigned int iter = 0;

		dual2 s = f(dual2(x, 1.0, 0.0));

		while(abs(s.Re()) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_CHEBYSHEV_ITER) {

			s = f(dual2(x, 1, 0));

			x = x - (s.Re() / s.Dual1())
						- square((s.Re() / s.Dual1()))
							* (s.Dual1() / (s.Dual2() * 2));
			iter++;
		}

		if(iter > OPTIMIZATION_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Approximate a root of a polynomial using Chebyshev's method.
	///
	/// @param p The polynomial to search the root of.
	/// @param guess The initial guess (defaults to 0).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline real root_chebyshev(polynomial<real> p, real guess = 0) {

		polynomial<> Dp = deriv(p);
		polynomial<> D2p = deriv(p);

		real x = guess;
		unsigned int iter = 0;

		while(abs(p(x)) > OPTIMIZATION_TOL && iter <= OPTIMIZATION_CHEBYSHEV_ITER) {
			x = x - (p(x) / Dp(x)) - square((p(x) / Dp(x))) * (Dp(x) / (D2p(x) * 2));
			iter++;
		}

		if(iter > OPTIMIZATION_CHEBYSHEV_ITER) {
			TH_MATH_ERROR("root_chebyshev", x, NO_ALGO_CONVERGENCE);
			return nan();
		}

		return x;
	}


	/// Find the roots of a function inside a given interval.
	///
	/// @param f The real function to search the root of.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @param steps The number of sub-intervals to check
	/// for alternating sign (optional).
	/// @return A vector of the roots of the function that were found.
	///
	/// @note If the number of roots inside the interval is completely unknown,
	/// using many more steps should be preferred, to ensure all roots are found.
	template<typename RealFunction>
	inline std::vector<real> roots(
		RealFunction f, real a, real b,
		real tolerance = OPTIMIZATION_TOL, real steps = 10) {

		if(steps == 0) {
			TH_MATH_ERROR("roots", steps, DIV_BY_ZERO);
			return {nan()};
		}

		// Find candidate intervals
		std::vector<vec2> intervals = find_root_intervals(f, a, b, steps);

		std::vector<real> res;
		res.reserve(intervals.size());

		// Iterate through all candidate intervals and refine the results
		for (unsigned int i = 0; i < intervals.size(); ++i) {

			// Check whether the extremes of the candidate intervals
			// happen to coincide with the roots
			if(abs(f(intervals[i][0])) <= MACH_EPSILON) {
				res.push_back(intervals[i][0]);
				continue;
			}

			if(abs(f(intervals[i][1])) <= MACH_EPSILON) {
				res.push_back(intervals[i][1]);
				continue;
			}

			// Approximate the roots using bisection inside each interval
			res.push_back(
				root_bisection(f, intervals[i][0], intervals[i][1], tolerance));
		}

		return res;
	}


	/// Find all the roots of a polynomial.
	/// An interval bound on the roots is found using Cauchy's theorem.
	///
	/// @param p The polynomial to search the roots of.
	/// @param steps The number of steps to use
	/// (defaults to twice the polynomial's order).
	/// @return A vector of the roots of the polynomial that were found.
	template<typename Field>
	inline std::vector<Field> roots(
		polynomial<Field> p, real tolerance = OPTIMIZATION_TOL,
		unsigned int steps = 0) {

		// Effective order of the polynomial
		const unsigned int n = p.find_order();
		p /= p.coeff[n];

		// Absolute value of the highest coefficient
		Field a_hi = abs(p.coeff[n]);
		Field a_sum = 0;

		// Sum the absolute values of the lesser coefficients
		for (unsigned int i = 0; i < n; ++i)
			a_sum += abs(p.coeff[i]);

		// The roots are bounded in absolute value by the maximum
		const Field M = max(a_hi, a_sum);

		// Back track from the bounds to the first sign inversion ?
		
		return roots(
			[p](real x) { return p(x); },
			-M, M, tolerance, steps ? steps : (2 * n));
	}

}

#endif
