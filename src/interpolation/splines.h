
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
				TH_MATH_ERROR("invlerp", t_new, MathError::OutOfDomain);
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
			TH_MATH_ERROR("slerp", P1_l ? P2_l : P1_l, MathError::ImpossibleOperation);
			return vec<real, N>(nan());
		}

		// Angle between P1 and P2 (from the dot product)
		real omega = acos((P1 * P2) / (P1_l * P2_l));
		real s = sin(omega);

		// Check that the sine of the angle is not zero
		if(s == 0) {
			TH_MATH_ERROR("slerp", s, MathError::DivByZero);
			return vec<real, N>(nan());
		}

		return (P1 * sin((1 - t) * omega) + P2 * sin(t * omega)) / s;
	}


	// Sigmoid-like interpolation


	/// Smoothstep interpolation
	inline real smoothstep(real x1, real x2, real interp) {

		if(x1 == x2) {
			TH_MATH_ERROR("smoothstep", x1, MathError::DivByZero);
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
			TH_MATH_ERROR("smootherstep", x1, MathError::DivByZero);
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
	inline vec<real, N> bezier_quadratic(
		const vec<real, N>& P0, const vec<real, N>& P1,
		const vec<real, N>& P2, real t) {

		return lerp(lerp(P0, P1, t), lerp(P1, P2, t), t);
	}


	/// Cubic Bezier curve
	template<unsigned int N>
	inline vec<real, N> bezier_cubic(
		const vec<real, N>& P0, const vec<real, N>& P1,
		const vec<real, N>& P2, vec<real, N> P3, real t) {

		const vec<real, N> A = lerp(P0, P1, t);
		const vec<real, N> B = lerp(P1, P2, t);
		const vec<real, N> C = lerp(P2, P3, t);

		const vec<real, N> D = lerp(A, B, t);
		const vec<real, N> E = lerp(B, C, t);

		return lerp(D, E, t);
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
	/// a set of data points using an optimized Thomas algorithm, with O(n) complexity.
	/// The x[i] values must be strictly increasing, an error is raised otherwise.
	///
	/// @param x The x[i] values of the data points.
	/// @param y The y[i] values of the data points.
	/// @return A vector of spline nodes representing the cubic spline interpolation.
	template<typename Dataset1, typename Dataset2>
	inline std::vector<spline_node> splines_cubic(
		const Dataset1& x, const Dataset2& y) {

		if(x.size() < 2) {
			TH_MATH_ERROR("splines_cubic", x.size(), MathError::InvalidArgument);
			return { spline_node(nan(), nan(), nan(), nan(), nan()) };
		}

		if(x.size() != y.size()) {
			TH_MATH_ERROR("splines_cubic", x.size(), MathError::InvalidArgument);
			return { spline_node(nan(), nan(), nan(), nan(), nan()) };
		}

		const size_t n_points = x.size();
		const size_t n_nodes = n_points - 1;

		if(n_points < 2) {
			TH_MATH_ERROR("splines_cubic", n_points, MathError::InvalidArgument);
			return { spline_node(nan(), nan(), nan(), nan(), nan()) };
		}

		// Check for strictly increasing x values
		for (size_t i = 0; i < n_nodes; ++i) {

			if(x[i + 1] <= x[i]) {
				TH_MATH_ERROR("splines_cubic", x[i + 1] - x[i], MathError::InvalidArgument);
				return { spline_node(nan(), nan(), nan(), nan(), nan()) };
			}
		}

		// Special case for linear interpolation
		if(n_points == 2) {
			const real slope = (y[1] - y[0]) / (x[1] - x[0]);
			return { spline_node(x[0], y[0], slope, 0, 0) };
		}

		// Second derivatives at each point
		std::vector<real> m (n_points);

		// Intervals, reused for diagonal ratios
		std::vector<real> h (n_nodes);
		
		// Compute differences h[i] = x[i+1] - x[i]
		for (size_t i = 0; i < n_nodes; ++i)
			h[i] = x[i + 1] - x[i];
		
		// Natural spline boundary conditions: m[0] = m[n] = 0
		m[0] = 0;
		m[n_nodes] = 0;
		
		// Solve tridiagonal system using Thomas algorithm
		// Am = b where m are the second derivatives,
		// using in-place forward elimination,
		// h array will store the modified diagonal ratios.
		
		// Compute RHS and setup for Thomas algorithm
		std::vector<real> rhs (n_points);
		rhs[0] = 0;
		rhs[n_nodes] = 0;
		
		// Build RHS: 6 * [(y[i+1] - y[i]) / h[i] - (y[i] - y[i-1]) / h[i-1]]
		for (size_t i = 1; i < n_nodes; ++i) {
			rhs[i] = 6.0 * (
				(y[i + 1] - y[i]) / h[i] -
				(y[i] - y[i - 1]) / h[i - 1]
			);
		}
		
		// Forward elimination

		// Tridiagonal matrix has form:
		// [ 2(h[0]+h[1])	h[1]			0,				...	]
		// [ h[1]			2(h[1]+h[2])	h[2]				]
		// [ 0				h[2]			2(h[2]+h[3])		]
		// [ 0				0				...				...	]
		
		// Diagonal elements
		std::vector<real> diag (n_points);

		// Boundary conditions
		diag[0] = 1.0;
		diag[n_nodes] = 1.0;
		
		for (size_t i = 1; i < n_nodes; ++i)
			diag[i] = 2.0 * (h[i - 1] + h[i]);
		
		// Forward sweep
		for (size_t i = 1; i < n_nodes; ++i) {
			const real ratio = h[i - 1] / diag[i - 1];
			diag[i] -= ratio * h[i - 1];
			rhs[i] -= ratio * rhs[i - 1];
		}
		
		// Back substitution
		for (int i = int(n_nodes) - 1; i >= 0; --i) {
			
			if(i < int(n_nodes) - 1)
				m[i] = (rhs[i] - h[i] * m[i + 1]) / diag[i];
			else
				m[i] = rhs[i] / diag[i];
		}

		// Compute spline coefficients from second derivatives
		// For the interval [x[i], x[i+1]], the spline is:
		// S(x) = a + b * (x - x[i]) + c * (x - x[i])^2 + d * (x - x[i])^3
		std::vector<spline_node> nodes (n_nodes);
		for (size_t i = 0; i < n_nodes; ++i) {
			
			// Coefficients in terms of second derivatives m and differences h
			const real a = y[i];
			const real b = (y[i + 1] - y[i]) / h[i] - h[i] * (2.0 * m[i] + m[i + 1]) / 6.0;
			const real c = m[i] * 0.5;
			const real d = (m[i + 1] - m[i]) / (6.0 * h[i]);
			
			nodes[i] = spline_node(x[i], a, b, c, d);
		}

		return nodes;
	}


	/// Compute the cubic splines interpolation of
	/// a set of data points using an optimized Thomas algorithm, with O(n) complexity.
	/// The x[i] values must be strictly increasing, an error is raised otherwise.
	///
	/// @param p The set of data points as a vector of coordinate pairs (x, y).
	/// @return A vector of spline nodes representing the cubic spline interpolation.
	template <typename DataPoints = std::vector<vec2>>
	inline std::vector<spline_node> splines_cubic(DataPoints p) {

		// Wrap a vector of vectors, providing access to a single coordinate as a vector
		struct accessor {

			const DataPoints& points;
			size_t j;

			accessor(const DataPoints& points, size_t j) : points(points), j(j) {}
			inline real operator[](size_t i) const { return points[i][j]; }
			inline size_t size() const { return points.size(); }
		};

		return splines_cubic(accessor(p, 0), accessor(p, 1));
	}


	/// @class spline
	/// A natural cubic spline interpolation class
	class spline {
	public:

		/// The computed nodes of the natural cubic spline
		/// interpolation over the points.
		std::vector<spline_node> nodes;


		/// Construct the natural cubic spline interpolation
		/// from a vector of coordinate pairs.
		///
		/// The x[i] values must be strictly increasing, an error is raised otherwise.
		template<typename DataPoints = std::vector<vec2>>
		spline(const DataPoints& p) {
			nodes = splines_cubic(p);
		}


		/// Construct the natural cubic spline interpolation
		/// from the sets of x[i] and y[i] data points.
		///
		/// The x[i] values must be strictly increasing, an error is raised otherwise.
		template<typename Dataset1, typename Dataset2>
		spline(const Dataset1& X, const Dataset2& Y) {
			nodes = splines_cubic(X, Y);
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
