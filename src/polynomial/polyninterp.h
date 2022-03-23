#ifndef UROBORO_POLYNINTERP_H
#define UROBORO_POLYNINTERP_H

#include <vector>
#include "./polynomial.h"
#include "../algebra/vec.h"


namespace uroboro {


	// Compute the Lagrange polynomial
	// interpolating a set of points
	template<typename T>
	inline polynomial<T> lagrange_polynomial(const std::vector<vec<2, T>>& points) {

		// Check that all x_i are different to prevent
		// division by zero
		for (int i = 0; i < points.size() - 1; ++i) {
			if(points[i].get(0) == points[i + 1].get(0)) {
				UMATH_ERROR("lagrange_polynomial", points[i].get(0), INVALID_ARGUMENT);
				return polynomial<T>({T(nan())});
			}
		}

		// Lagrange polynomial to construct
		polynomial<T> L = {0};

		for (int j = 0; j < points.size(); ++j) {
			
			// The Lagrange polynomial is a linear
			// combination of all l_j
			polynomial<T> l_j = {1};

			for (int m = 0; m < points.size(); ++m) {
				
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


	// Compute the interpolating polynomial of a real function
	// on an equidistant point sample.
	// <f> is the function to interpolate
	// <a> and <b> are the extremes of the interval (a < b)
	// <order> is the order of the resulting polynomial
	inline polynomial<real> interpolate_grid(real_function f, real a, real b, unsigned int order) {

		std::vector<vec2> points;
		points.resize(order + 1);

		// Sample <order + 1> equidistant points
		for (int i = 0; i < order + 1; ++i) {
			real x = (b - a) / real(order) * real(i);
			points[i] = {x, f(x)};
		}

		return lagrange_polynomial(points);
	}

}

#endif
