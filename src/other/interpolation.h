#ifndef UROBORO_INTERP_H
#define UROBORO_INTERP_H

#include "../constants.h"

namespace uroboro {


	// Linear interpolation
	inline real lerp(real x1, real x2, real interp) {
		return (x1 + interp * (x2 - x1));
	}


	// Linear interpolation
	template<unsigned int N>
	inline vec<N> lerp(vec<N> P1, vec<N> P2, real interp) {
		return (P1 + (P2 - P1) * interp);
	}


	// Inverse linear interpolation
	inline real invlerp(real x1, real x2, real value) {
		return (value - x1) / (x2 - x1);
	}


	// Inverse linear interpolation
	template<unsigned int N>
	inline vec<N> invlerp(vec<N> P1, vec<N> P2, real value) {

		real t = (value - P1.get(0)) / (P2.get(0) - P1.get(0));

		// Check that all computed t_i are the same
		for (int i = 1; i < N; ++i) {
			if(t != (value - P1.get(i)) / (P2.get(i) - P1.get(i))) {
				UMATH_ERROR("invlerp",
					(value - P1.get(i)) / (P2.get(i) - P1.get(i)),
					OUT_OF_DOMAIN);
				return nan();
			}
		}

		return t;
	}


	// Remap a value from one range to another
	inline real remap(real iFrom, real iTo, real oFrom, real oTo, real value) {
		return lerp(oFrom, oTo, invlerp(iFrom, iTo, value));
	}


	// Remap a vector value from one range to another
	template<unsigned int N>
	inline vec<N> remap(vec<N> iFrom, vec<N> iTo, vec<N> oFrom, vec<N> oTo, real value) {
		return lerp(oFrom, oTo, invlerp(iFrom, iTo, value));
	}


	// Sigmoid-like interpolation


	// Smoothstep interpolation
	inline real smoothstep(real x1, real x2, real interp) {

		// Clamp x between 0 and 1
		const real x = clamp((interp - x1) / (x2 - x1), 0.0, 1.0);

		// 3x^2 - 2x^3
		return x * x * (3 - 2 * x);
	}


	// Smootherstep interpolation
	inline real smootherstep(real x1, real x2, real interp) {

		// Clamp x between 0 and 1
		const real x = clamp((interp - x1) / (x2 - x1), 0.0, 1.0);

		// 6x^5 - 15x^4 + 10x^3
		return x * x * x * (x * (x * 6 - 15) + 10);
	}


	// Bezier curves


	// Quadratic Bezier curve
	template<unsigned int N>
	inline vec<N> quadratic_bezier(vec<N> P0, vec<N> P1, vec<N> P2, real t) {
		return lerp(lerp(P0, P1, t), lerp(P1, P2, t), t);
	}


	// Cubic Bezier curve
	template<unsigned int N>
	inline vec<N> cubic_bezier(vec<N> P0, vec<N> P1, vec<N> P2, vec<N> P3, real t) {

		vec<N> A = lerp(P0, P1, t);
		vec<N> B = lerp(P1, P2, t);
		vec<N> C = lerp(P2, P3, t);

		vec<N> D = lerp(A, B, t);
		vec<N> E = lerp(B, C, t);

		return lerp(D, E, t);
	}
	
}

#endif
