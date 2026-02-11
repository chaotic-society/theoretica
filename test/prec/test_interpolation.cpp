
///
/// @file test_interpolation.cpp Test cases for interpolation functions and splines
///

#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


// Runge's function
real runge(real x) {
	return 1.0 / (1.0 + 25.0 * x * x);
}

// Cubic polynomial
real cubic(real x) {
	return x * x * x - 2 * x + 1;
}


// Generate a random vector with gaussian elements
vec2 rand_vec2(random::random_source& rnd) {
	return vec2({
		rnd.gaussian(0, 1E+08), rnd.gaussian(0, 1E+08)
	});
}

// Distance function for polynomial comparison
long double distance_polyn(const polynomial<real>& p1, const polynomial<real>& p2) {
	const polynomial<real> d = p1 - p2;
	real r = -inf();

	for (size_t i = 0; i < d.size(); ++i)
		r = max(r, th::abs(d[i]));

	return r;
}

// Distance function for vector comparison
template<unsigned int N>
long double distance_vec(const vec<real, N>& v1, const vec<real, N>& v2) {
	return (v1 - v2).norm();
}


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("interpolation", argc, argv);
	ctx.output->settings.outputFiles = { "test/prec/prec_interpolation.csv" };
	random::random_source rnd = ctx.random->get_rnd();

	auto polyn_opt = prec::equation_options<polynomial<real>>(
		1E-08,
		distance_polyn
	);

	auto vec2_opt = prec::equation_options<vec2>(
		1E-08,
		distance_vec<2>
	);


	// polynomial.h


	// Test lagrange_polynomial with simple points
	{
		std::vector<vec2> points = {{0, 1}, {1, 2}, {2, 5}};
		polynomial<real> p = lagrange_polynomial(points);
		
		// Polynomial should pass through all points
		ctx.equals("lagrange_polynomial(point[0])", p(0.0), 1.0);
		ctx.equals("lagrange_polynomial(point[1])", p(1.0), 2.0);
		ctx.equals("lagrange_polynomial(point[2])", p(2.0), 5.0);
	}

	// Test lagrange_polynomial with quadratic polynomial
	{
		// y = x^2 + x + 1
		std::vector<vec2> points = {{0, 1}, {1, 3}, {2, 7}};
		polynomial<real> p = lagrange_polynomial(points);
		polynomial<real> expected = {1.0, 1.0, 1.0};  // 1 + x + x^2
		
		ctx.equals("lagrange_polynomial(quadratic)", p, expected, polyn_opt);
	}

	// Test lagrange_polynomial with cubic polynomial
	{
		// y = x^3 - 2x + 1
		std::vector<vec2> points = {{-1, 2}, {0, 1}, {1, 0}, {2, 5}};
		polynomial<real> p = lagrange_polynomial(points);
		
		ctx.equals("lagrange_polynomial(cubic[0])", p(-1.0), 2.0);
		ctx.equals("lagrange_polynomial(cubic[1])", p(0.0), 1.0);
		ctx.equals("lagrange_polynomial(cubic[2])", p(1.0), 0.0);
		ctx.equals("lagrange_polynomial(cubic[3])", p(2.0), 5.0);
	}

	// Test lagrange_polynomial interpolates exactly for polynomials
	{
		std::vector<vec2> points = {{0, 1}, {1, 2}, {2, 9}, {3, 28}};
		polynomial<real> p = lagrange_polynomial(points);
		
		// Should match x^3 + 1 at intermediate points
		ctx.equals("lagrange_polynomial(interpolation)", p(1.5), 1.5 * 1.5 * 1.5 + 1.0, 1E-06);
	}

	// Test chebyshev_nodes
	{
		auto nodes = chebyshev_nodes(0.0, 1.0, 5);
		
		ctx.equals("chebyshev_nodes (size)", (int) nodes.size(), 5);
		
		// All nodes should be in [0, 1]
		for (size_t i = 0; i < nodes.size(); ++i) {
			
			ctx.equals("chebyshev_nodes (range)", (nodes[i] >= 0.0 && nodes[i] <= 1.0), true);
		}
	}

	// Test chebyshev_nodes symmetry
	{
		auto nodes = chebyshev_nodes(-1.0, 1.0, 10);
		
		// Nodes should be symmetric around 0
		for (size_t i = 0; i < nodes.size() / 2; ++i) {

			ctx.equals("chebyshev_nodes (symmetry)", 
				std::abs(nodes[i] + nodes[nodes.size() - 1 - i]),
				0.0, 1E-10
			);
		}
	}

	// Test interpolate_grid on polynomial (should be exact)
	{
		polynomial<real> p = interpolate_grid(cubic, -2.0, 2.0, 3);
		
		// Should exactly match the cubic polynomial
		ctx.equals("interpolate_grid(cubic[0])", p(0.0), cubic(0.0), 1E-10);
		ctx.equals("interpolate_grid(cubic[1])", p(1.0), cubic(1.0), 1E-10);
		ctx.equals("interpolate_grid(cubic[2])", p(-1.5), cubic(-1.5), 1E-10);
	}

	// Test interpolate_grid on sine
	{
		polynomial<real> p = interpolate_grid(CAST_LAMBDA(std::sin, real), 0.0, PI / 2, 5);
		
		// Check at interpolation points
		for (int i = 0; i <= 5; ++i) {
			real x = (PI / 2) / 5.0 * i;
			ctx.equals("interpolate_grid(sin)", p(x), std::sin(x), 1E-06);
		}
	}

	// Test interpolate_chebyshev on polynomial (should be exact)
	{
		polynomial<real> p = interpolate_chebyshev(cubic, -2.0, 2.0, 3);
		
		ctx.equals("interpolate_chebyshev(cubic)", p(0.5), cubic(0.5), 1E-08);
	}


	// splines.h


	// Test lerp
	{
		ctx.equals("lerp(0.0)", lerp(0.0, 10.0, 0.0), 0.0);
		ctx.equals("lerp(0.5)", lerp(0.0, 10.0, 0.5), 5.0);
		ctx.equals("lerp(1.0)", lerp(0.0, 10.0, 1.0), 10.0);
		ctx.equals("lerp(0.25)", lerp(5.0, 15.0, 0.25), 7.5);
	}

	// Test lerp vector
	{
		vec2 p1 = rand_vec2(rnd);
		vec2 p2 = rand_vec2(rnd);
		
		ctx.equals("lerp(vec, 0.0)", lerp(p1, p2, 0.0), p1, vec2_opt);
		ctx.equals("lerp(vec, 1.0)", lerp(p1, p2, 1.0), p2, vec2_opt);
		ctx.equals("lerp(vec, 0.5)", lerp(p1, p2, 0.5), (p2 + p1) * 0.5, vec2_opt);
	}

	// Test invlerp
	{
		ctx.equals("invlerp(start)", invlerp(0.0, 10.0, 0.0), 0.0);
		ctx.equals("invlerp(mid)", invlerp(0.0, 10.0, 5.0), 0.5);
		ctx.equals("invlerp(end)", invlerp(0.0, 10.0, 10.0), 1.0);
		ctx.equals("invlerp(quarter)", invlerp(5.0, 15.0, 7.5), 0.25);
	}

	// Test remap
	{
		ctx.equals("remap(identity)", remap(0.0, 10.0, 0.0, 10.0, 5.0), 5.0);
		ctx.equals("remap(scale)", remap(0.0, 10.0, 0.0, 100.0, 5.0), 50.0);
		ctx.equals("remap(shift)", remap(0.0, 10.0, 10.0, 20.0, 5.0), 15.0);
		ctx.equals("remap(general)", remap(0.0, 1.0, -1.0, 1.0, 0.5), 0.0);
	}

	// Test lerp and invlerp are inverses
	{
		real x1 = 5.0, x2 = 15.0, val = 8.0;
		real t = invlerp(x1, x2, val);
		ctx.equals("invlerp(lerp)", lerp(x1, x2, t), val);
	}

	// Test nlerp
	{
		vec<real, 3> p1 = {1.0, 0.0, 0.0};
		vec<real, 3> p2 = {0.0, 1.0, 0.0};
		vec<real, 3> result = nlerp(p1, p2, 0.5);
		
		// Result should be normalized
		ctx.equals("nlerp (norm)", result.norm(), 1.0, 1E-10);
	}

	// Test slerp maintains length
	{
		vec<real, 3> p1 = {1.0, 0.0, 0.0};
		vec<real, 3> p2 = {0.0, 1.0, 0.0};
		vec<real, 3> result = slerp(p1, p2, 0.5);
		
		ctx.equals("slerp (norm)", result.norm(), 1.0, 1E-08);
	}

	// Test slerp endpoints
	{
		vec<real, 3> p1 = {1.0, 0.0, 0.0};
		vec<real, 3> p2 = {0.0, 1.0, 0.0};
		
		auto slerp_opt = prec::equation_options<vec<real, 3>>(
			1E-08,
			[](const vec<real, 3>& v1, const vec<real, 3>& v2) { return (v1 - v2).norm(); }
		);
		
		ctx.equals("slerp(t = 0)", slerp(p1, p2, 0.0), p1, slerp_opt);
		ctx.equals("slerp(t = 1)", slerp(p1, p2, 1.0), p2, slerp_opt);
	}

	// Test smoothstep boundary conditions
	{
		ctx.equals("smoothstep(0.0)", smoothstep(0.0, 1.0, 0.0), 0.0);
		ctx.equals("smoothstep(1.0)", smoothstep(0.0, 1.0, 1.0), 1.0);
		ctx.equals("smoothstep(0.5)", smoothstep(0.0, 1.0, 0.5), 0.5);
	}

	// Test smoothstep is smooth (derivative is 0 at endpoints)
	{
		real h = 1E-08;
		real deriv_0 = (smoothstep(0.0, 1.0, h) - smoothstep(0.0, 1.0, 0.0)) / h;
		real deriv_1 = (smoothstep(0.0, 1.0, 1.0) - smoothstep(0.0, 1.0, 1.0 - h)) / h;
		
		ctx.equals("smoothstep(0 +- h)", deriv_0, 0.0, 1E-06);
		ctx.equals("smoothstep(1 +- h)", deriv_1, 0.0, 1E-06);
	}

	// Test smoothstep clamping
	{
		ctx.equals("smoothstep(clamp < 0)", smoothstep(0.0, 1.0, -0.5), 0.0);
		ctx.equals("smoothstep(clamp > 1)", smoothstep(0.0, 1.0, 1.5), 1.0);
	}

	// Test smootherstep boundary conditions
	{
		ctx.equals("smootherstep(0.0)", smootherstep(0.0, 1.0, 0.0), 0.0);
		ctx.equals("smootherstep(1.0)", smootherstep(0.0, 1.0, 1.0), 1.0);
		ctx.equals("smootherstep(0.5)", smootherstep(0.0, 1.0, 0.5), 0.5);
	}

	// Test smootherstep clamping
	{
		ctx.equals("smootherstep(clamp<0)", smootherstep(0.0, 1.0, -0.5), 0.0);
		ctx.equals("smootherstep(clamp>1)", smootherstep(0.0, 1.0, 1.5), 1.0);
	}

	// Test quadratic_bezier endpoints
	{
		vec2 p0 = {0.0, 0.0};
		vec2 p1 = {0.5, 1.0};
		vec2 p2 = {1.0, 0.0};
		
		ctx.equals("quadratic_bezier(t = 0)", quadratic_bezier(p0, p1, p2, 0.0), p0, vec2_opt);
		ctx.equals("quadratic_bezier(t = 1)", quadratic_bezier(p0, p1, p2, 1.0), p2, vec2_opt);
	}

	// Test quadratic_bezier midpoint
	{
		vec2 p0 = {0.0, 0.0};
		vec2 p1 = {1.0, 2.0};
		vec2 p2 = {2.0, 0.0};
		vec2 mid = quadratic_bezier(p0, p1, p2, 0.5);
		
		vec2 expected = {1.0, 1.0};
		ctx.equals("quadratic_bezier(t = 0.5)", mid, expected, vec2_opt);
	}

	// Test cubic_bezier endpoints
	{
		vec2 p0 = {0.0, 0.0};
		vec2 p1 = {0.33, 1.0};
		vec2 p2 = {0.66, 1.0};
		vec2 p3 = {1.0, 0.0};
		
		ctx.equals("cubic_bezier(t = 0)", cubic_bezier(p0, p1, p2, p3, 0.0), p0, vec2_opt);
		ctx.equals("cubic_bezier(t = 1)", cubic_bezier(p0, p1, p2, p3, 1.0), p3, vec2_opt);
	}

	// Test generic bezier with 2 points (should be lerp)
	{
		std::vector<vec2> points = {{0.0, 0.0}, {1.0, 1.0}};
		vec2 result = bezier(points, 0.5);
		vec2 expected = lerp(points[0], points[1], 0.5);
		
		ctx.equals("bezier(2 points)", result, expected, vec2_opt);
	}

	// Test generic bezier endpoints
	{
		std::vector<vec2> points = {{0.0, 0.0}, {0.5, 1.0}, {1.0, 0.0}};
		
		ctx.equals("bezier(t = 0)", bezier(points, 0.0), points[0], vec2_opt);
		ctx.equals("bezier(t = 1)", bezier(points, 1.0), points.back(), vec2_opt);
	}

	// Test cubic_splines interpolates points
	{
		std::vector<vec2> points = {{0, 0}, {1, 1}, {2, 4}, {3, 9}};
		std::vector<theoretica::spline_node> nodes = cubic_splines(points);

		// Check interpolation at data points
		for (size_t i = 0; i < points.size() - 1; ++i) {
			ctx.equals("cubic_splines(points)", nodes[i](points[i][0]), points[i][1], 1E-10);
		}
	}

	// Test cubic_splines with separate X and Y
	{
		std::vector<real> x = {0, 1, 2, 3};
		std::vector<real> y = {0, 1, 4, 9};
		auto nodes = cubic_splines(x, y);
		
		// Check interpolation
		for (size_t i = 0; i < x.size() - 1; ++i) {
			ctx.equals("cubic_splines(X, Y)", nodes[i](x[i]), y[i], 1E-10);
		}
	}

	// Test spline class construction
	{
		std::vector<vec2> points = {{0, 0}, {1, 1}, {2, 0}, {3, -1}};
		spline s(points);
		
		// Check interpolation at data points
		for (const auto& p : points) {
			ctx.equals("spline(points)", s(p[0]), p[1], 1E-10);
		}
	}

	// Test spline with separate datasets
	{
		std::vector<real> x = {0, 1, 2, 3, 4};
		std::vector<real> y = {0, 1, 0, -1, 0};
		spline s(x, y);
		
		for (size_t i = 0; i < x.size(); ++i) {
			ctx.equals("spline(X, Y)", s(x[i]), y[i], 1E-10);
		}
	}

	// Test spline continuity
	{
		std::vector<vec2> points = {{0, 0}, {1, 1}, {2, 0}};
		spline s(points);
		
		// Evaluate just before and after knot
		real h = 1E-06;
		real before = s(1.0 - h);
		real after = s(1.0 + h);
		
		ctx.equals("spline (continuity)", before, after, 1E-04);
	}

	// Test spline derivative
	{
		std::vector<vec2> points = {{0, 0}, {1, 1}, {2, 4}, {3, 9}};
		spline s(points);
		
		// Derivative should be continuous
		real h = 1E-06;
		real deriv_before = s.deriv(1.5 - h);
		real deriv_after = s.deriv(1.5 + h);
		
		ctx.equals("spline.deriv() (continuity)", deriv_before, deriv_after, 1E-04);
	}

	// Test spline on known function
	{
		// Interpolate sin(x) on [0, pi]
		std::vector<vec2> points;
		for (int i = 0; i <= 10; ++i) {
			real x = PI / 10.0 * i;
			points.push_back({x, std::sin(x)});
		}
		
		spline s(points);
		
		// Check accuracy at intermediate points
		auto opt = prec::estimate_options<real, real>(
			prec::interval(0.0, PI),
			prec::estimator::quadrature1D(),
			1E-02, 1000
		);
		
		ctx.estimate("spline(sin)",
			[&s](real x) { return s(x); },
			[](real x) { return std::sin(x); },
			opt
		);
	}
	

	// Estimator tests for interpolation accuracy


	// Test Lagrange grid interpolation on smooth functions
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(-0.8, 0.8),
			prec::estimator::quadrature1D(),
			1E-01, ctx.settings.defaultIterations
		);
		
		polynomial<real> p = interpolate_grid([](real x) {
			return std::cos(x);
		}, -1.0, 1.0, 12);
		
		ctx.estimate("interpolate_grid(cos, 12)",
			[&p](real x) { return p(x); },
			CAST_LAMBDA(std::cos, real),
			opt
		);
	}

	// Test Lagrange interpolation with Chebyshev nodes on smooth functions
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(-0.8, 0.8),
			prec::estimator::quadrature1D(),
			1E-01, ctx.settings.defaultIterations
		);
		
		polynomial<real> p = interpolate_chebyshev([](real x) {
			return std::cos(x);
		}, -1.0, 1.0, 12);
		
		ctx.estimate("interpolate_chebyshev(cos, 12)",
			[&p](real x) { return p(x); },
			CAST_LAMBDA(std::cos, real),
			opt
		);
	}

	// Test Chebyshev interpolation is better than grid for Runge
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(-0.8, 0.8),
			prec::estimator::quadrature1D(),
			1E-1, ctx.settings.defaultIterations
		);
		
		polynomial<real> p = interpolate_chebyshev(runge, -1.0, 1.0, 16);
		
		ctx.estimate("interpolate_chebyshev(runge, 12)",
			[&p](real x) { return p(x); },
			runge,
			opt
		);
	}
}
