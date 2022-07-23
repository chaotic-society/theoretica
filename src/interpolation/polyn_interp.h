
///
/// @file polyn_interp.h Polynomial interpolation of real functions
///

#ifndef THEORETICA_POLYN_INTERP_H
#define THEORETICA_POLYN_INTERP_H

#include <vector>
#include "../core/real_analysis.h"
#include "../polynomial/polynomial.h"
#include "../algebra/vec.h"
#include "../core/function.h"


namespace theoretica {


	/// Compute the Lagrange polynomial
	/// interpolating a set of points
	/// @param points The set of n points to interpolate
	/// @return A polynomial of (n - 1) degree interpolating the points
	template<typename T = real>
	inline polynomial<T> lagrange_polynomial(const std::vector<vec<2, T>>& points) {

		if(!points.size()) {
			TH_MATH_ERROR("lagrange_polynomial", points.size(), INVALID_ARGUMENT);
			return polynomial<T>({T(nan())});
		}

		// Check that all x_i are different to prevent
		// division by zero
		for (unsigned int i = 0; i < points.size() - 1; ++i) {
			if(points[i].get(0) == points[i + 1].get(0)) {
				TH_MATH_ERROR("lagrange_polynomial", points[i].get(0), INVALID_ARGUMENT);
				return polynomial<T>({T(nan())});
			}
		}

		// Lagrange polynomial to construct
		polynomial<T> L = {0};

		for (unsigned int j = 0; j < points.size(); ++j) {
			
			// The Lagrange polynomial is a linear
			// combination of all l_j
			polynomial<T> l_j = {1};

			for (unsigned int m = 0; m < points.size(); ++m) {
				
				if(m == j)
					continue;

				// l_j = product(x - x_m / x_j - x_m)
				l_j *= polynomial<T>({-points[m].get(0), 1});
				l_j /= points[j].get(0) - points[m].get(0);
			}

			// L = sum(y_j * l_j)
			l_j *= points[j].get(1);
			L += l_j;
		}

		return L;
	}


	/// Compute the n Chebyshev nodes on a given interval
	/// @param a The lower bound of the interval
	/// @param b The upper bound of the interval
	/// @param n The number of points to evaluate
	inline vec_buff chebyshev_nodes(real a, real b, unsigned int n) {

		vec_buff nodes;
		nodes.resize(n);

		real m = (b + a) / 2.0;
		real c = (b - a) / 2.0;

		for (unsigned int i = 1; i < n + 1; ++i)
			nodes[i - 1] = m + c * cos(real(2 * i - 1) / real(2 * n) * PI);

		return nodes;
	}


	/// Compute the interpolating polynomial of a real function
	/// on an equidistant point sample.
	/// @param f The function to interpolate
	/// @param a Lower bound of the interval
	/// @param b Upper bound of the interval 
	/// @param order Order of the resulting polynomial
	/// @return A polynomial of (n - 1) degree interpolating the function
	inline polynomial<real> interpolate_grid(real_function f, real a, real b, unsigned int order) {

		std::vector<vec2> points;
		points.resize(order + 1);

		// Sample <order + 1> equidistant points
		for (unsigned int i = 0; i < order + 1; ++i) {
			real x = (b - a) / real(order) * real(i);
			points[i] = {x, f(x)};
		}

		return lagrange_polynomial(points);
	}


	/// Compute the interpolating polynomial of a real function
	/// using Chebyshev nodes as sampling points
	/// @param f The function to interpolate
	/// @param a Lower bound of the interval
	/// @param b Upper bound of the interval 
	/// @param order Order of the resulting polynomial
	/// @return A polynomial of (n - 1) degree interpolating the function through
	/// the Chebyshev nodes.
	///
	/// \see chebyshev_nodes
	/// \see lagrange_polynomial
	inline polynomial<real> interpolate_chebyshev(real_function f, real a, real b, unsigned int order) {

		std::vector<vec2> points;
		points.resize(order + 1);

		vec_buff nodes = chebyshev_nodes(a, b, order + 1);

		// Sample <order + 1> equidistant points
		for (unsigned int i = 0; i < order + 1; ++i)
			points[i] = {nodes[i], f(nodes[i])};

		return lagrange_polynomial(points);
	}

}

#endif
