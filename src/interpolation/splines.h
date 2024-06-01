
///
/// @file splines.h Spline interpolation
///

#ifndef THEORETICA_SPLINES_H
#define THEORETICA_SPLINES_H

#include "../core/constants.h"
#include "../algebra/algebra_types.h"


namespace theoretica {


	/// Linear interpolation
	inline real lerp(real x1, real x2, real interp) {
		return (x1 + interp * (x2 - x1));
	}


	/// Linear interpolation
	template<unsigned int N>
	inline vec<real, N> lerp(
		const vec<real, N>& P1, const vec<real, N>& P2, real interp) {

		return (P1 + (P2 - P1) * interp);
	}


	/// Inverse linear interpolation
	inline real invlerp(real x1, real x2, real value) {
		return (value - x1) / (x2 - x1);
	}


	/// Inverse linear interpolation
	template<unsigned int N>
	inline vec<real, N> invlerp(
		const vec<real, N>& P1, const vec<real, N>& P2, real value) {

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
	inline vec<real, N> remap(
		const vec<real, N>& iFrom, const vec<real, N>& iTo,
		const vec<real, N>& oFrom, const vec<real, N>& oTo, real value) {
		return lerp(oFrom, oTo, invlerp(iFrom, iTo, value));
	}


	/// Normalized linear interpolation
	template<unsigned int N>
	inline vec<real, N> nlerp(
		const vec<real, N>& P1, const vec<real, N>& P2, real interp) {
		
		return (P1 + (P2 - P1) * interp).normalized();
	}


	/// Spherical interpolation
	template<unsigned int N>
	inline vec<real, N> slerp(
		const vec<real, N>& P1, const vec<real, N>& P2, real t) {

		// Compute (only once) the length
		// of the input vectors
		const real P1_l = P1.norm();
		const real P2_l = P2.norm();

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
	inline vec<real, N> quadratic_bezier(
		const vec<real, N>& P0, const vec<real, N>& P1,
		const vec<real, N>& P2, real t) {

		return lerp(lerp(P0, P1, t), lerp(P1, P2, t), t);
	}


	/// Cubic Bezier curve
	template<unsigned int N>
	inline vec<real, N> cubic_bezier(
		const vec<real, N>& P0, const vec<real, N>& P1,
		const vec<real, N>& P2, vec<real, N> P3, real t) {

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
	inline vec<real, N> bezier(const std::vector<vec<real, N>>& points, real t) {
		
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


	/// A cubic splines node for a given x interval
	struct spline_node {

		/// Upper extreme of the interpolation interval \f$x_i\f$
		real x;

		/// Coefficients of the interpolating cubic spline
		/// \f$a + b x + c x^2 + d x^3\f$
		real a, b, c, d;

		/// Default constructor
		spline_node() : x(0), a(0), b(0), c(0), d(0) {}

		/// Construct from \f$x_i\f$ and polynomial coefficients
		spline_node(real x, real a, real b, real c, real d)
			: x(x), a(a), b(b), c(c), d(d) {}


		/// Evaluate the interpolating cubic spline
		/// (no check on the input value is performed!)
		inline real operator()(real X) const {
			const real h = X - x;
			return a + h * (b + h * (c + h * d));
		}


		/// Evaluate the derivative of the interpolating cubic spline
		/// (no check on the input value is performed)
		inline real deriv(real X) const {
			const real h = X - x;
			return b + h * (c * 2 + h * d * 3);
		}
	};


	/// Compute the cubic splines interpolation of
	/// a set of data points.
	/// The X values should be strictly increasing.
	template<typename DataPoints = std::vector<vec2>>
	inline std::vector<spline_node> cubic_splines(DataPoints p) {

		if(p.size() < 2) {
			TH_MATH_ERROR("cubic_splines", p.size(), INVALID_ARGUMENT);
			return {spline_node(nan(), nan(), nan(), nan(), nan())};
		}

		const unsigned int n = p.size() - 1;

		std::vector<real> dx(n);
		std::vector<real> delta(n);

		delta[0] = 0;

		std::vector<real> alpha(n + 1);
		std::vector<real> beta(n + 1);
		std::vector<real> gamma(n + 1);

		alpha[0] = 1;
		beta[0] = 0;
		gamma[0] = 0;

		std::vector<real> b(n);
		std::vector<real> c(n);
		std::vector<real> d(n);

		for (unsigned int i = 0; i < n; ++i)
			dx[i] = p[i + 1][0] - p[i][0];

		for (unsigned int i = 1; i < n; ++i)
			delta[i] = 3 * (
				((p[i + 1][1] - p[i][1]) / dx[i])
					- (p[i][1] - p[i - 1][1]) / dx[i - 1]);

		for (unsigned int i = 1; i < n; ++i) {
			alpha[i] = 2 * (p[i + 1][0] - p[i - 1][0]) - dx[i - 1] * beta[i - 1];
			beta[i] = dx[i] / alpha[i];
			gamma[i] = (delta[i] - dx[i] * gamma[i - 1]) / alpha[i];
		}

		// Apply boundary conditions
		alpha[n] = 1;
		gamma[n] = 0;
		c[n] = 0;

		// Solve the associated tridiagonal system
		// using back-substitution
		for (int i = n - 1; i >= 0; --i) {
			
			c[i] = gamma[i] - beta[i] * c[i + 1];
			b[i] = (p[i + 1][1] - p[i][1]) / dx[i]
				- dx[i] * (c[i + 1] + 2 * c[i]) / 3.0;
			d[i] = (c[i + 1] - c[i]) / (3.0 * dx[i]);
		}

		std::vector<spline_node> nodes(n);

		for (unsigned int i = 0; i < n; ++i)
			nodes[i] = spline_node(p[i][0], p[i][1], b[i], c[i], d[i]);

		return nodes;
	}


	/// Compute the cubic splines interpolation of
	/// the sets of X and Y data points.
	/// The X values should be strictly increasing.
	template<typename Dataset1, typename Dataset2>
	inline std::vector<spline_node> cubic_splines(
		const Dataset1& x, const Dataset2& y) {

		if(x.size() < 2) {
			TH_MATH_ERROR("cubic_splines", x.size(), INVALID_ARGUMENT);
			return {spline_node(nan(), nan(), nan(), nan(), nan())};
		}

		if(x.size() != y.size()) {
			TH_MATH_ERROR("cubic_splines", x.size(), INVALID_ARGUMENT);
			return {spline_node(nan(), nan(), nan(), nan(), nan())};
		}

		const unsigned int n = x.size() - 1;

		std::vector<real> dx(n);
		std::vector<real> delta(n);

		delta[0] = 0;

		std::vector<real> alpha(n + 1);
		std::vector<real> beta(n + 1);
		std::vector<real> gamma(n + 1);

		alpha[0] = 1;
		beta[0] = 0;
		gamma[0] = 0;

		std::vector<real> b(n);
		std::vector<real> c(n);
		std::vector<real> d(n);

		for (unsigned int i = 0; i < n; ++i)
			dx[i] = x[i + 1] - x[i];

		for (unsigned int i = 1; i < n; ++i)
			delta[i] = 3 * (
				((y[i + 1] - y[i]) / dx[i])
					- (y[i] - y[i - 1]) / dx[i - 1]);

		for (unsigned int i = 1; i < n; ++i) {
			alpha[i] = 2 * (x[i + 1] - x[i - 1]) - dx[i - 1] * beta[i - 1];
			beta[i] = dx[i] / alpha[i];
			gamma[i] = (delta[i] - dx[i] * gamma[i - 1]) / alpha[i];
		}

		// Apply boundary conditions
		alpha[n] = 1;
		gamma[n] = 0;
		c[n] = 0;

		// Solve the associated tridiagonal system
		// using back-substitution
		for (int i = n - 1; i >= 0; --i) {
			
			c[i] = gamma[i] - beta[i] * c[i + 1];
			b[i] = (y[i + 1] - y[i]) / dx[i]
				- dx[i] * (c[i + 1] + 2 * c[i]) / 3.0;
			d[i] = (c[i + 1] - c[i]) / (3.0 * dx[i]);
		}

		std::vector<spline_node> nodes(n);

		for (unsigned int i = 0; i < n; ++i)
			nodes[i] = spline_node(x[i], y[i], b[i], c[i], d[i]);

		return nodes;
	}


	/// @class spline
	/// A natural cubic spline interpolation class
	class spline {
	public:

		/// The computed nodes of the natural cubic spline
		/// interpolation over the points.
		std::vector<spline_node> nodes;


		/// Construct the natural cubic spline interpolation
		/// from a set of X and y data points.
		template<typename DataPoints = std::vector<vec2>>
		spline(const DataPoints& p) {
			nodes = cubic_splines(p);
		}


		/// Construct the natural cubic spline interpolation
		/// from the sets of X and Y data points.
		template<typename Dataset1, typename Dataset2>
		spline(const Dataset1& X, const Dataset2& Y) {
			nodes = cubic_splines(X, Y);
		}


		/// Copy from another spline
		inline spline& operator=(const spline& other) {
			nodes = other.nodes;
			return *this;
		}


		/// Evaluate the natural cubic spline interpolation
		/// at a given point.
		inline real operator()(real x) const {
			
			for (int i = int(nodes.size()) - 1; i > 0; --i)
				if(x >= nodes[i].x)
					return nodes[i](x);

			// Extrapolation for x < x_0
			return nodes[0](x);
		}


		/// Evaluate the derivative of the natural cubic
		/// spline interpolation at a given point.
		inline real deriv(real x) const {
			
			for (unsigned int i = nodes.size() - 1; i > 0; --i)
				if(x >= nodes[i].x)
					return nodes[i].deriv(x);

			// Extrapolation for x < x_0
			return nodes[0].deriv(x);
		}


		/// Get an iterator to the first spline element.
		inline auto begin() {
			return nodes.begin();
		}


		/// Get an iterator to one plus the last spline element.
		inline auto end() {
			return nodes.end();
		}
	};

	
}

#endif
