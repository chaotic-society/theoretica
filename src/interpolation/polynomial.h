
///
/// @file polynomial.h Polynomial interpolation of real functions
///

#ifndef THEORETICA_INTERP_POLYNOMIAL_H
#define THEORETICA_INTERP_POLYNOMIAL_H

#include <vector>
#include "../core/real_analysis.h"
#include "../polynomial/polynomial.h"
#include "../algebra/algebra_types.h"
#include "../core/function.h"


namespace theoretica {


	/// Compute the Lagrange interpolating polynomial of a set of points.
	/// This function allows for interpolation of any type that supports the
	/// operations used in the Lagrange formula (multiplication and division
	/// by scalars, addition), including vector y points.
	///
	/// @param x The x coordinates of the points
	/// @param y The y coordinates of the points
	/// @return A polynomial of (n - 1) degree interpolating the points
	/// @tparam Vector1 A vector type for the x coordinates
	/// @tparam Vector2 A vector type for the y coordinates
	template <
		typename Vector1, typename Vector2,
		enable_vector<Vector1> = true,
		enable_vector<Vector2> = true,
		typename Type = vector_element_t<Vector2>
	>
	inline polynomial<Type> lagrange(const Vector1& x, const Vector2& y) {

		if(!x.size()) {
			TH_MATH_ERROR("lagrange", x.size(), MathError::InvalidArgument);
			return make_error<polynomial<Type>>(1);
		}

		if(x.size() != y.size()) {
			TH_MATH_ERROR("lagrange", y.size(), MathError::InvalidArgument);
			return make_error<polynomial<Type>>(1);
		}

		// Check that all x_i are different to prevent division by zero
		for (unsigned int i = 0; i < x.size() - 1; ++i) {

			if(abs(x[i] - x[i + 1]) < MACH_EPSILON) {
				TH_MATH_ERROR("lagrange", x[i], MathError::InvalidArgument);
				return make_error<polynomial<Type>>(1);
			}
		}

		// Lagrange polynomial to construct
		polynomial<Type> L = {0};

		for (unsigned int j = 0; j < x.size(); ++j) {
			
			// The Lagrange polynomial is a linear combination of all l_j
			polynomial<real> l_j = {1.0};

			for (unsigned int m = 0; m < x.size(); ++m) {
				
				if(m == j)
					continue;

				// l_j = product(x - x_m / x_j - x_m)
				l_j *= polynomial<real>({-x[m], 1.0}) / (x[j] - x[m]);
			}

			// L = sum(y_j * l_j)
			polynomial<Type> term;
			term.resize(l_j.size());

			for (size_t i = 0; i < l_j.size(); ++i)
				term[i] = y[j] * l_j[i];

			L += term;
		}

		return L;
	}


	/// Compute the n Chebyshev nodes on a given interval. These nodes
	/// are the roots of the Chebyshev polynomial of degree n, and are
	/// optimal for polynomial interpolation, minimizing Runge's phenomenon.
	///
	/// @param a The lower bound of the interval
	/// @param b The upper bound of the interval
	/// @param n The number of points to evaluate
	/// @return A vector of the n Chebyshev nodes in the interval [a, b]
	template<typename Vector = vec<real>>
	inline Vector chebyshev_nodes(real a, real b, unsigned int n) {

		Vector nodes;
		nodes.resize(n);

		real m = (b + a) / 2.0;
		real c = (b - a) / 2.0;

		for (unsigned int i = 1; i < n + 1; ++i)
			nodes[i - 1] = m + c * cos(real(2 * i - 1) / real(2 * n) * PI);

		return nodes;
	}


	/// Compute the interpolating polynomial of a real function
	/// on an equidistant grid.
	///
	/// @param f The function to interpolate
	/// @param a Lower bound of the interval
	/// @param b Upper bound of the interval 
	/// @param order Order of the resulting polynomial
	/// @return A polynomial of (n - 1) degree interpolating the function
	template <
		typename RealFunction,
		typename Type = return_type_t<RealFunction>
	>
	inline polynomial<Type> interpolate_grid(RealFunction f, real a, real b, unsigned int order) {

		vec<real> x (order + 1);
		vec<Type> y (order + 1);

		// Sample <order + 1> equidistant points
		for (unsigned int i = 0; i < order + 1; ++i) {
			x[i] = (b - a) / real(order) * real(i) + a;
			y[i] = f(x[i]);
		}

		return lagrange(x, y);
	}


	/// Compute the interpolating polynomial of a real function
	/// using Chebyshev nodes as sampling points.
	///
	/// @param f The function to interpolate
	/// @param a Lower bound of the interval
	/// @param b Upper bound of the interval 
	/// @param order Order of the resulting polynomial
	/// @return A polynomial of (n - 1) degree interpolating the function through
	/// the Chebyshev nodes.
	///
	/// \see chebyshev_nodes
	/// \see lagrange
	template <
		typename RealFunction,
		typename Type = return_type_t<RealFunction>
	>
	inline polynomial<Type> interpolate_chebyshev(RealFunction f, real a, real b, unsigned int order) {

		vec<real> x = chebyshev_nodes(a, b, order + 1);
		vec<Type> y (order + 1);

		// Sample <order + 1> equidistant points
		for (unsigned int i = 0; i < order + 1; ++i)
			y[i] = f(x[i]);

		return lagrange(x, y);
	}


	/// Use the best available method to compute the interpolating polynomial of a real function.
	///
	/// @param f The function to interpolate
	/// @param a Lower bound of the interval
	/// @param b Upper bound of the interval
	/// @param order Order of the resulting polynomial
	/// @return A polynomial of (n - 1) degree interpolating the function
	template <
		typename RealFunction,
		typename Type = return_type_t<RealFunction>
	>
	inline polynomial<Type> interpolate(RealFunction f, real a, real b, unsigned int order) {
		return interpolate_chebyshev(f, a, b, order);
	}

}

#endif
