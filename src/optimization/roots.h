
///
/// @file roots.h Root approximation of real functions
///

#ifndef THEORETICA_ROOTS_H
#define THEORETICA_ROOTS_H

#include "../core/function.h"
#include "../calculus/deriv.h"
#include "../autodiff/dual.h"
#include "../autodiff/dual2.h"
#include "../algebra/vec.h"
#include "../complex/complex.h"
#include "../core/iter_result.h"

#include <iostream>

namespace theoretica {


	/// Find candidate intervals for root finding by evaluating a function
	/// at equidistant points inside an interval [a, b] and checking its sign.
	///
	/// @param f A function of real variable
	/// @param a The lower extreme of the interval
	/// @param b The upper extreme of the interval
	/// @param steps The number of sub-intervals to check (defaults to 10)
	template<typename RealFunction, typename Bracket = vec2>
	inline std::vector<Bracket> find_root_intervals(
		RealFunction f, real a, real b, unsigned int steps = 10) {

		std::vector<Bracket> res;
		const real dx = (b - a) / (real) steps;

		for (unsigned int i = 0; i < steps; ++i) {
			
			const real x1 = a + i * dx;
			const real x2 = a + (i + 1) * dx;

			if(f(x1) * f(x2) <= 0)
				res.push_back({x1, x2});
		}

		return res;
	}


	/// Find the root of a univariate real function using bisection
	/// inside a compact interval [a, b] where \f$f(a) * f(b) < 0\f$.
	///
	/// @param f The univariate real function.
	/// @param a The lower extreme of the interval.
	/// @param b The upper extreme of the interval.
	/// @param tol The minimum half-length of the bracketing interval to stop
	/// the algorithm, so that \f$|x_r - x_l| \leq 2\epsilon\f$.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline iter_result<real> root_bisect(
		RealFunction f, real a, real b,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_BISECTION_ITER) {

		if(a > b) {
			TH_MATH_ERROR("root_bisect", a, MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		if(f(a) * f(b) >= 0) {
			TH_MATH_ERROR("root_bisect", f(a) * f(b), MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		real x_avg = 0.0;
		real x_min = a;
		real x_max = b;

		unsigned int iter = 0;

		while((x_max - x_min) > (2 * tol) && iter <= max_iter) {

			x_avg = (x_max + x_min) / 2.0;

			if(f(x_avg) * f(x_min) > 0)
				x_min = x_avg;
			else
				x_max = x_avg;

			++iter;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_bisect", x_avg, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(b - a) / 2.0);
		}

		return iter_result<real>(x_avg, iter, abs(b - a) / 2.0);
	}


	/// Find a root of a univariate real function using the ITP
	/// (Interpolate-Truncate-Project) method, by bracketing the zero
	/// inside a compact interval [a, b] where \f$f(a) * f(b) < 0\f$.
	/// The \f$k_2\f$ parameter is chosen to be 2, avoiding expensive
	/// operations while retaining good convergence. This method is
	/// the best choice when the function is not smooth and is expensive
	/// to compute.
	///
	/// @param f The univariate real function.
	/// @param a The lower extreme of the interval.
	/// @param b The upper extreme of the interval.
	/// @param tol The minimum half-length of the bracketing interval to stop
	/// the algorithm, so that \f$|x_r - x_l| \leq 2\epsilon\f$.
	/// @param n0 A hyper-parameter of the algorithm, must be zero or greater. Bigger values
	/// give more importance to the regula falsi estimate.
	/// @param k1 A hyper-parameter of the algorithm, influences the truncation step
	/// (defaults to \f$0.2 / (b - a)\f$, following the authors' results).
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline real root_itp(
		RealFunction f, real a, real b, real tol = OPTIMIZATION_TOL,
		unsigned int n0 = 1, real k1 = 0.0) {

		// Default value for k1
		if (k1 == 0.0)
			k1 = 0.2 / (b - a);

		if(a > b) {
			TH_MATH_ERROR("root_itp", a, MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		real y_a = f(a);
		real y_b = f(b);

		if(y_a * y_b >= 0) {
			TH_MATH_ERROR("root_itp", y_a * y_b, MathError::InvalidArgument);
			return iter_result<real>(ConvergenceStatus::InvalidInput);
		}

		// Monotonicity of the function
		const int monotone = (y_a < y_b) ? 1 : -1;

		real x_t, x_new;

		const long int n_half = th::floor(th::log2((b - a) / tol));
		const long int n_max = n_half + n0;

		real eps = tol * pow(2.0, n_max);
		long int iter = 0;

		while((b - a) > (2 * tol) && iter <= n_max) {


			// Interpolation
			const real x_f = (a * y_b - b * y_a) / (y_b - y_a);
			const real x_half = 0.5 * (a + b);


			// Truncation
			const int sigma = sgn(x_half - x_f);
			const real delta = k1 * square(b - a);

			if (delta <= abs(x_half - x_f))
				x_t = x_f + sigma * delta;
			else
				x_t = x_half;


			// Projection
			const real r = eps - (b - a) / 2.0;

			if (abs(x_t - x_half) <= r)
				x_new = x_t;
			else
				x_new = x_half - sigma * r;


			// Update
			const real y_new = f(x_new);

			if (monotone * y_new > 0) {
				b = x_new;
				y_b = y_new;
			} else if (monotone * y_new < 0) {
				a = x_new;
				y_a = y_new;
			} else {
				return x_new;
			}

			eps *= 0.5;
			iter++;
		}

		if(abs(b - a) > 2 * tol) {
			TH_MATH_ERROR("root_itp", (a + b) / 2.0, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(b - a) / 2.0);
		}

		return iter_result<real>((a + b) / 2.0, iter, abs(b - a) / 2.0);
	}


	/// Find a root of a univariate real function using Newton's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline iter_result<real> root_newton(
		RealFunction f, RealFunction Df, real guess = 0.0,
		real tol = OPTIMIZATION_TOL, unsigned int max_iter = OPTIMIZATION_NEWTON_ITER) {

		real x = guess;
		real f_x = inf();
		real Df_x = 0;
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			Df_x = Df(x);

			// Check for division by zero
			if (abs(Df_x) < MACH_EPSILON) {
				TH_MATH_ERROR("root_newton", Df_x, MathError::DivByZero);
				return iter_result<real>(ConvergenceStatus::Diverged, iter, f_x);
			}

			f_x = f(x);
			x = x - (f_x / Df_x);
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_newton", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find a root of a univariate real function using Newton's method,
	/// computing the derivative using automatic differentiation.
	///
	/// @param f The real function to search the root of,
	/// with dual argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline iter_result<real> root_newton(
		dual(*f)(dual), real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_NEWTON_ITER) {

		real x = guess;
		dual s = dual(inf(), 0.0);
		unsigned int iter = 0;

		while(abs(s.Re()) > tol && iter <= max_iter) {

			// Compute the function and its derivative at the same time
			s = f(dual(x, 1.0));

			x = x - (s.Re() / s.Dual());
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_newton", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(s.Re()));
		}

		return iter_result<real>(x, iter, abs(s.Re()));
	}


	/// Find a root of a complex function using Newton's method.
	///
	/// @param f The complex function to search the root of.
	/// @param df The derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(z_n)|^2 \leq \epsilon^2\f$.
	/// @param max_iter The maximum number of iterations before stopping
	/// the algorithm (defaults to OPTIMIZATION_NEWTON_ITER).
	/// @return The coordinate of the root of the function,
	/// or a complex NaN if the algorithm did not converge.
	template <
		typename Type = real,
		typename ComplexFunction = std::function<complex<Type>(complex<Type>)>
	>
	inline iter_result<complex<Type>> root_newton(
		ComplexFunction f, ComplexFunction Df, complex<Type> guess,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_NEWTON_ITER) {

		complex<Type> z = guess;
		complex<Type> f_z = Type(inf());
		unsigned int iter = 0;

		while(f_z * f_z.conjugate() > tol * tol && iter <= max_iter) {

			f_z = f(z);
			z = z - (f_z / df(z));
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_newton", z.Re(), MathError::NoConvergence);
			return iter_result<complex<Type>>(ConvergenceStatus::MaxIterations, iter, f_z.norm());
		}

		return iter_result<complex<Type>>(z, iter, f_z.norm());
	}


	/// Find a root of a univariate real function using Halley's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline iter_result<real> root_halley(
		RealFunction f, RealFunction Df,
		RealFunction D2f, real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_HALLEY_ITER) {

		real x = guess;
		real f_x = inf();
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			f_x = f(x);
			const real df_x = Df(x);

			x = x - (2 * f_x * df_x) / (2 * square(df_x) - f_x * D2f(x));
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_halley", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find a root of a univariate real function using Halley's method,
	/// leveraging automatic differentiation to compute the first and second
	/// derivatives of the function.
	///
	/// @param f The real function to search the root of,
	/// with dual2 argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| < \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	inline iter_result<real> root_halley(
		dual2(*f)(dual2), real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_HALLEY_ITER) {

		real x = guess;
		dual2 s = dual2(inf(), 0.0, 0.0);
		unsigned int iter = 0;

		while(abs(s.Re()) > tol && iter <= max_iter) {

			// Compute the function value and the first and
			// second derivatives at the same time.
			s = f(dual2(x, 1, 0));

			const real f_x = s.Re();
			const real df_x = s.Dual1();
			const real d2f_x = s.Dual2();

			x = x - (2 * f_x * df_x) / (2 * square(df_x) - f_x * d2f_x);
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_halley", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(s.Re()));
		}

		return iter_result<real>(x, iter, abs(s.Re()));
	}


	/// Find a root of a univariate real function using Steffensen's method.
	///
	/// @param f The real function to search the root of.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	template<typename RealFunction>
	inline iter_result<real> root_steffensen(
		RealFunction f, real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_STEFFENSEN_ITER) {

		real x = guess;
		real f_x = inf();
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			const real f_x = f(x);
			const real g_x = (f(x + f_x) / f_x) - 1;

			x = x - (f_x / g_x);
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_steffensen", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find a root of a univariate real function using Chebyshev's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param D2f The second derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	///
	/// Chebyshev's method can be derived by expanding the inverse of the function
	/// around the zero and truncating the series. This method is particularly suited
	/// when the derivatives of the function are easy to compute, especially
	/// when using automatic differentiation.
	template<typename RealFunction>
	inline real root_chebyshev(
		RealFunction f, RealFunction Df,
		RealFunction D2f, real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_CHEBYSHEV_ITER) {

		real x = guess;
		real f_x = inf();
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			f_x = f(x);
			const real df_x = Df(x);
			const real u = f_x / df_x;

			x = x - u - square(u) * D2f(x) / (2.0 * df_x);

			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_chebyshev", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find a root of a univariate real function using Chebyshev's method,
	/// by computing the first and second derivatives using automatic differentiation.
	///
	/// @param f The real function to search the root of,
	/// with dual2 argument and return value.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	///
	/// Chebyshev's method can be derived by expanding the inverse of the function
	/// around the zero and truncating the series. This method is particularly suited
	/// when the derivatives of the function are easy to compute, especially
	/// when using automatic differentiation. For example, the exponential needs to be
	/// computed only once to evaluate the function and its derivatives.
	inline iter_result<real> root_chebyshev(
		dual2(*f)(dual2), real guess = 0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_CHEBYSHEV_ITER) {

		real x = guess;
		dual2 s = dual2(inf(), 0.0, 0.0);
		unsigned int iter = 0;

		while(abs(s.Re()) > tol && iter <= max_iter) {

			s = f(dual2(x, 1.0, 0.0));

			const real f_x = s.Re();
			const real df_x = s.Dual1();
			const real u = f_x / df_x;

			x = x - u - square(u) * s.Dual2() / (2.0 * df_x);
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_chebyshev", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(s.Re()));
		}

		return iter_result<real>(x, iter, abs(s.Re()));
	}


	/// Find a root of a univariate real function using Ostrowski's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	///
	/// Ostrowski's method is a 4-th order method using 2 function evaluations
	/// and 1 function evaluation. It combines a Newton step with a
	/// corrective coefficient. This method does not have an overload using
	/// automatic differentiation as it would hinder the performance benefit.
	template<typename RealFunction>
	inline iter_result<real> root_ostrowski(
		RealFunction f, RealFunction Df, real guess = 0.0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_OSTROWSKI_ITER) {

		real x = guess;
		real f_x = inf();
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			f_x = f(x);
			const real df_x = Df(x);
			const real u = f_x / df_x;
			const real f_xu = f(x - u);

			x = x - u - (f_xu / df_x) * (f_x / (f_x - 2 * f_xu));

			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_ostrowski", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find a root of a univariate real function using Jarrat's method.
	///
	/// @param f The real function to search the root of.
	/// @param Df The first derivative of the function.
	/// @param guess The initial guess (defaults to 0).
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|f(x_n)| \leq \epsilon\f$.
	/// @param max_iter The maximum number of iterations to perform, after which
	/// the algorithm is assumed to not have converged.
	/// @return The coordinate of the root of the function,
	/// or NaN if the algorithm did not converge.
	///
	/// Jarrat's method is a 4-th order method which is particularly suited when
	/// the derivative of the function is less computationally expensive than
	/// the function itself, like in the case of integrals. This method does not
	/// have an overload using automatic differentiation as it would hinder the
	/// performance benefit.
	template<typename RealFunction>
	inline iter_result<real> root_jarrat(
		RealFunction f, RealFunction Df, real guess = 0.0,
		real tol = OPTIMIZATION_TOL,
		unsigned int max_iter = OPTIMIZATION_JARRAT_ITER) {

		real x = guess;
		real f_x = inf();
		unsigned int iter = 0;

		while(abs(f_x) > tol && iter <= max_iter) {

			f_x = f(x);
			const real df_x = Df(x);
			const real u = f_x / df_x;
			const real f_xu = Df(x - 2.0 * u / 3.0);

			x = x - 0.5 * u + f_x / (df_x - 3 * f_xu);
			iter++;
		}

		if(iter > max_iter) {
			TH_MATH_ERROR("root_jarrat", x, MathError::NoConvergence);
			return iter_result<real>(ConvergenceStatus::MaxIterations, iter, abs(f_x));
		}

		return iter_result<real>(x, iter, abs(f_x));
	}


	/// Find the roots of a univariate real function inside a given interval,
	/// by first searching for candidate intervals and then applying bracketing methods.
	///
	/// @param f The real function to search the root of.
	/// @param a The lower extreme of the search interval.
	/// @param b The upper extreme of the search interval.
	/// @param tol The \f$\epsilon\f$ tolerance value to stop the algorithm
	/// when \f$|x_r - x_l| \leq 2\epsilon\f$.
	/// @param steps The number of sub-intervals to check
	/// for alternating sign (defaults to 10).
	/// @return A vector of the roots of the function that were found.
	///
	/// @note If the number of roots inside the interval is completely unknown,
	/// using many more steps should be preferred, to ensure all roots are found.
	template<typename RealFunction, typename Vector = vec<real>>
	inline Vector roots(
		RealFunction f, real a, real b,
		real tol = OPTIMIZATION_TOL, real steps = 10) {

		if(steps == 0) {
			TH_MATH_ERROR("roots", steps, MathError::DivByZero);
			return {nan()};
		}

		// Find candidate intervals
		auto intervals = find_root_intervals(f, a, b, steps);

		Vector res;
		res.resize(intervals.size());
		size_t idx = 0;

		// Iterate through all candidate intervals and refine the results
		for (unsigned int i = 0; i < intervals.size(); ++i) {

			// Check whether the extremes of the candidate intervals
			// happen to coincide with the roots
			if(abs(f(intervals[i][0])) <= MACH_EPSILON) {
				res[idx] = intervals[i][0];
				idx++;
				continue;
			}

			if(abs(f(intervals[i][1])) <= MACH_EPSILON) {
				res[idx] = intervals[i][1];
				idx++;
				continue;
			}

			// Approximate the roots using bisection inside each interval
			real r = root_bisect(f, intervals[i][0], intervals[i][1], tol);

			if (!is_nan(r)) {
				res[idx] = r;
				idx++;
			}
		}

		// Remove extra elements
		//res.resize(idx);
		return res;
	}


	/// Find all the roots of a polynomial.
	/// An interval bound on the roots is found using Cauchy's theorem.
	///
	/// For polynomials with very unevenly distributed roots, the algorithm may
	/// fail to find all roots, in this case the steps parameter should be increased
	/// to ensure all roots are found.
	///
	/// @param p The polynomial to search the roots of.
	/// @param steps The number of steps to use (defaults to 10 times the polynomial's order).
	/// @return A vector of the roots of the polynomial that were found.
	template<typename Field, typename Vector = vec<Field>>
	inline Vector roots_polynomial(
		const polynomial<Field>& p, real tolerance = OPTIMIZATION_TOL,
		unsigned int steps = 0) {
		
		// Real polynomial degree
		const unsigned int p_order = p.find_order();
		if (p_order == 0) {
			TH_MATH_ERROR("roots_polynomial", p.coeff[0], MathError::InvalidArgument);
			return {nan()};
		}

		// Normalize the polynomial by the leading coefficient to apply Cauchy's bound
		polynomial<Field> p_norm = p / p.coeff[p_order];
		Field a_max = 0;

		// Find the maximum absolute coefficient
		for (unsigned int i = 0; i < p_order - 1; ++i)
			a_max = max(a_max, abs(p_norm.coeff[i]));

		// The roots are bounded in absolute value by the maximum
		const Field M = 1 + a_max;

		return roots(p_norm, -M, M, tolerance, steps > 0 ? steps : 10 * p_order);
	}

}

#endif
