
///
/// @file spline_interp.h Spline interpolation
///

#ifndef THEORETICA_INTERP_H
#define THEORETICA_INTERP_H

#include "../core/constants.h"


namespace theoretica {


	/// Linear interpolation
	inline real lerp(real x1, real x2, real interp) {
		return (x1 + interp * (x2 - x1));
	}


	/// Linear interpolation
	template<unsigned int N>
	inline vec<real, N> lerp(vec<real, N> P1, vec<real, N> P2, real interp) {
		return (P1 + (P2 - P1) * interp);
	}


	/// Inverse linear interpolation
	inline real invlerp(real x1, real x2, real value) {
		return (value - x1) / (x2 - x1);
	}


	/// Inverse linear interpolation
	template<unsigned int N>
	inline vec<real, N> invlerp(vec<real, N> P1, vec<real, N> P2, real value) {

		real t = invlerp(P1.get(0), P2.get(0), value);

		// Check that all computed t_i are the same
		for (int i = 1; i < N; ++i) {

			real t_new = invlerp(P1.get(i), P2.get(i), value);

			if(t != t_new) {
				TH_MATH_ERROR("invlerp", t_new, OUT_OF_DOMAIN);
				return nan();
			}
		}

		return t;
	}


	/// Remap a value from one range to another
	inline real remap(real iFrom, real iTo, real oFrom, real oTo, real value) {
		return lerp(oFrom, oTo, invlerp(iFrom, iTo, value));
	}


	/// Remap a vector value from one range to another
	template<unsigned int N>
	inline vec<real, N> remap(vec<real, N> iFrom, vec<real, N> iTo, vec<real, N> oFrom, vec<real, N> oTo, real value) {
		return lerp(oFrom, oTo, invlerp(iFrom, iTo, value));
	}


	/// Normalized linear interpolation
	template<unsigned int N>
	inline vec<real, N> nlerp(vec<real, N> P1, vec<real, N> P2, real interp) {
		return (P1 + (P2 - P1) * interp).normalized();
	}


	/// Spherical interpolation
	template<unsigned int N>
	inline vec<real, N> slerp(vec<real, N> P1, vec<real, N> P2, real t) {

		// Compute (only once) the length
		// of the input vectors
		real P1_l = P1.norm();
		real P2_l = P2.norm();

		// Check whether one of the vectors is null,
		// which would make the computation impossible
		if(P1_l == 0 || P2_l == 0) {
			TH_MATH_ERROR("slerp", P1_l ? P2_l : P1_l, IMPOSSIBLE_OPERATION);
			return vec<real, N>(nan());
		}

		// Angle between P1 and P2 (from the dot product)
		real omega = acos((P1 * P2) / (P1_l * P2_l));
		real s = sin(omega);

		// Check that the sine of the angle is not zero
		if(s == 0) {
			TH_MATH_ERROR("slerp", s, DIV_BY_ZERO);
			return vec<real, N>(nan());
		}

		return (P1 * sin((1 - t) * omega) + P2 * sin(t * omega)) / s;
	}


	// Sigmoid-like interpolation


	/// Smoothstep interpolation
	inline real smoothstep(real x1, real x2, real interp) {

		if(x1 == x2) {
			TH_MATH_ERROR("smoothstep", x1, DIV_BY_ZERO);
			return nan();
		}

		// Clamp x between 0 and 1
		const real x = clamp((interp - x1) / (x2 - x1), 0.0, 1.0);

		// 3x^2 - 2x^3
		return x * x * (3 - 2 * x);
	}


	/// Smootherstep interpolation
	inline real smootherstep(real x1, real x2, real interp) {

		if(x1 == x2) {
			TH_MATH_ERROR("smootherstep", x1, DIV_BY_ZERO);
			return nan();
		}

		// Clamp x between 0 and 1
		const real x = clamp((interp - x1) / (x2 - x1), 0.0, 1.0);

		// 6x^5 - 15x^4 + 10x^3
		return x * x * x * (x * (x * 6 - 15) + 10);
	}


	// Bezier curves


	/// Quadratic Bezier curve
	template<unsigned int N>
	inline vec<real, N> quadratic_bezier(vec<real, N> P0, vec<real, N> P1, vec<real, N> P2, real t) {
		return lerp(lerp(P0, P1, t), lerp(P1, P2, t), t);
	}


	/// Cubic Bezier curve
	template<unsigned int N>
	inline vec<real, N> cubic_bezier(vec<real, N> P0, vec<real, N> P1, vec<real, N> P2, vec<real, N> P3, real t) {

		vec<real, N> A = lerp(P0, P1, t);
		vec<real, N> B = lerp(P1, P2, t);
		vec<real, N> C = lerp(P2, P3, t);

		vec<real, N> D = lerp(A, B, t);
		vec<real, N> E = lerp(B, C, t);

		return lerp(D, E, t);
	}


	/// Generic Bezier curve in N dimensions
	/// @param points The control points
	/// @param t The curve parameter between 0 and 1
	///
	/// The generic Bezier curve is computed by
	/// successive linear interpolations.
	/// For cubic and quadratic Bezier curves the
	/// related functions should be preferred.
	/// \see quadratic_bezier
	/// \see cubic_bezier
	template<unsigned int N>
	inline vec<real, N> bezier(std::vector<vec<real, N>> points, real t) {
		
		if(points.size() < 2) {
			TH_MATH_ERROR("bezier", points.size(), INVALID_ARGUMENT);
			return vec<real, N>(nan());
		}

		if(t < 0 || t > 1) {
			TH_MATH_ERROR("bezier", t, INVALID_ARGUMENT);
			return vec<real, N>(nan());
		}

		for (int index = points.size(); index > 1; --index) {

			for (int i = 0; i < index - 1; ++i)
				points[i] = lerp(points[i], points[i + 1], t);
		}

		return points[0];
	}
	
}

#endif
