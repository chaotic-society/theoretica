
#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

using namespace chebyshev;
using namespace theoretica;
using namespace autodiff;


// Simple univariate functions for testing

real poly2(real x) { return x * x - 3 * x + 2; }
real dpoly2(real x) { return 2 * x - 3; }

real poly3(real x) { return x * x * x + 2 * x * x - x + 5; }
real dpoly3(real x) { return 3 * x * x + 4 * x - 1; }
real d2poly3(real x) { return 6 * x + 4; }


// Harmonic functions (Laplacian = 0)

// h1: ln(x^2 + y^2) - harmonic away from origin
dual2 h1(vec<dual2> v) {
	return ln(v[0] * v[0] + v[1] * v[1]);
}

// h2: e^x * sin(y) - harmonic 
dual2 h2(vec<dual2> v) {
	return exp(v[0]) * sin(v[1]);
}

// h3: e^(-y) * cos(x) - harmonic
dual2 h3(vec<dual2> v) {
	return exp(-v[1]) * cos(v[0]);
}


// Multidual test functions

// Distance function for matrices
template<typename Matrix = mat<real>>
auto mat_distance = [](const Matrix& m1, const Matrix& m2) -> long double {

	if (m1.rows() != m2.rows() || m1.cols() != m2.cols())
		return inf();
	
	long double d = 0.0;
	for (size_t i = 0; i < m1.rows(); ++i)
		for (size_t j = 0; j < m1.cols(); ++j)
			d = th::max(d, th::abs(m1(i, j) - m2(i, j)));

	return d;
};


int main(int argc, char const *argv[]) {

	auto ctx = prec::make_context("autodiff", argc, argv);
	ctx.settings.outputFiles = { "test/prec/prec_autodiff.csv" };
	auto rnd = ctx.random->get_rnd();

	const real MAX = 1E+03;

	// multidual.h

	// Test gradient of a linear function
	{
		auto f = [](dvec v) { return 2 * v[0] - 3 * v[1] + v[2]; };
		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		vec<real> expected = {2.0, -3.0, 1.0};
		auto grad = gradient(f, point);

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance, prec::distance::euclidean<vec<real>>
		);

		ctx.equals("gradient(linear, [1,2,3])",
			grad,
			expected,
			opt
		);
	}

	// Test gradient of quadratic form ||v||^2
	{
		auto f = [](dvec v) { return v * v; };

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		vec<real> expected = {2.0 * point[0], 2.0 * point[1]};
		auto grad = gradient(f, point);

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance, prec::distance::euclidean<vec<real>>
		);

		ctx.equals("gradient(v^2, [3,4])",
			grad,
			expected,
			opt
		);
	}

	// Test gradient of product v[0]*v[1]
	{
		auto f = [](dvec v) { return v[0] * v[1]; };
		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		vec<real> expected = {point[1], point[0]};
		auto grad = gradient(f, point);

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance, prec::distance::euclidean<vec<real>>
		);

		ctx.equals("gradient(v0*v1, [2,5])",
			grad,
			expected,
			opt
		);
	}

	// Test gradient of 3D quadratic
	{
		auto f = [](dvec v) { return v[0] * v[0] + v[1] * v[1] + v[2] * v[2]; };
		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		vec<real> expected = {2.0 * point[0], 2.0 * point[1], 2.0 * point[2]};
		auto grad = gradient(f, point);

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance, prec::distance::euclidean<vec<real>>
		);

		ctx.equals(
			"gradient(x2+y2+z2, [1,2,3])",
			grad, expected, opt
		);
	}

	// Test gradient linearity
	{
		auto f = [](dvec v) { return v[0] * v[0]; };
		auto g = [](dvec v) { return v[1] * v[1]; };
		auto sum = [](dvec v) { return v[0] * v[0] + v[1] * v[1]; };

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto grad_f = gradient(f, point);
		auto grad_g = gradient(g, point);
		auto grad_sum = gradient(sum, point);

		vec<real> expected = grad_f + grad_g;
		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance, prec::distance::euclidean<vec<real>>
		);

		ctx.equals(
			"d(f+g) = df+dg",
			grad_sum, expected, opt
		);
	}

	// Test Laplacian of harmonic function h1: ln(x^2 + y^2) ≈ 0
	{
		vec<real> point = {
			rnd.uniform(0.5, MAX),
			rnd.uniform(0.5, MAX)
		};

		real lap = laplacian(h1, point);

		// ln(r^2) is harmonic away from origin
		ctx.equals("laplacian(ln(x2+y2), [1,1])",
			lap, 0.0
		);
	}

	// Test Laplacian of harmonic function h2: e^x * sin(y) ≈ 0
	{
		vec<real> point = {
			rnd.uniform(-100.0, 10.0),
			rnd.uniform(-th::PI, th::PI)
		};

		real lap = laplacian(h2, point);

		// e^x * sin(y) is harmonic: d2/dx2 + d2/dy2 = e^x*sin(y) - e^x*sin(y) = 0
		ctx.equals("laplacian(exp(x)*sin(y), [0.5,0.5])",
			lap,
			0.0
		);
	}

	// Test Laplacian of harmonic function h3: e^(-y) * cos(x)
	{
		vec<real> point = {
			rnd.uniform(-th::PI, th::PI),
			rnd.uniform(-100.0, 10.0)
		};

		real lap = laplacian(h3, point);

		// e^(-y)*cos(x) is harmonic: d2/dx2 + d2/dy2 = -e^(-y)*cos(x) + e^(-y)*cos(x) = 0
		ctx.equals("laplacian(exp(-y)*cos(x), [0.3,0.7])",
			lap,
			0.0
		);
	}

	// Test Laplacian of function: x^2 + 2*y^2
	{
		auto f = [](vec<dual2> v) { 
			return v[0] * v[0] + 2 * v[1] * v[1]; 
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real lap = laplacian(f, point);

		// laplacian(x^2 + 2*y^2) = 2 + 4 = 6
		ctx.equals("laplacian(x2+2y2, [1,2])", lap, 6.0);
	}

	// Test divergence of identity field in 2D: div([x,y]) = 1 + 1 = 2
	{
		auto f = [](dvec v) -> dvec {
			return { v[0], v[1] };
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real div = divergence(f, point);

		ctx.equals("divergence(identity_2D)", div, 2.0);
	}

	// Test divergence of scaled field [2x,2y,2z]
	{
		auto f = [](dvec v) {
			return 2.0 * v;
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real div = divergence(f, point);

		ctx.equals("divergence(2v)", div, 2.0 * point.size());
	}

	// Test divergence-free field: F = [-y, x, 0] -> div(F) = 0
	{
		auto f = [](dvec v) -> dvec {
			
			return {
				-v[1], v[0], dreal(0.0)
			};
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real div = divergence(f, point);

		ctx.equals("divergence([-y,x,0])", div, 0.0);
	}

	// Test divergence of gradient of x^2 + y^2: laplacian(x^2 + y^2) = 4
	{
		auto f = [](vec<dual2> v) { 
			return v[0] * v[0] + v[1] * v[1];
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real lap = laplacian(f, point);

		// laplacian(x^2+y^2) = 2 + 2 = 4
		ctx.equals("laplacian(x2+y2)", lap, 4.0);
	}

	// Test divergence of gradient of x^2 + 2*y^2 + 3*z^2: laplacian(...) = 2 + 4 + 6 = 12
	{
		auto f = [](vec<dual2> v) { 
			return v[0] * v[0] + 2.0 * v[1] * v[1] + 3.0 * v[2] * v[2];
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		real lap = laplacian(f, point);

		// laplacian(x^2 + 2y^2 + 3z^2) = 2 + 4 + 6 = 12
		ctx.equals("laplacian(x^2+2y^2+3z^2)", lap, 12.0);
	}

	// Test divergence (Laplacian) of harmonic function: laplacian(ln(x^2+y^2)) = 0
	{
		auto f = [](vec<dual2> v) {
			return ln(v[0] * v[0] + v[1] * v[1]);
		};

		vec<real> point = {
			rnd.uniform(0.5, MAX),
			rnd.uniform(0.5, MAX)
		};

		real lap = laplacian(f, point);

		// ln(x^2+y^2) is harmonic
		ctx.equals("laplacian(ln(r2))", lap, 0.0);
	}

	// Test divergence of harmonic function: laplacian(e^x * sin(y)) = 0
	{
		auto f = [](vec<dual2> v) {
			return exp(v[0]) * sin(v[1]);
		};

		vec<real> point = {
			rnd.uniform(-100.0, 10.0),
			rnd.uniform(-th::PI, th::PI)
		};
		
		real lap = laplacian(f, point);

		// e^x * sin(y) is harmonic
		ctx.equals("laplacian(exp(x)*sin(y))",
			lap,
			0.0
		);
	}

	// Test Jacobian of identity field: J = I
	{
		auto f = [](dvec v) { return v; };

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto J = jacobian(f, point);

		mat<real> expected = algebra::identity<mat<real>>(2, 2);

		auto opt = prec::equation_options<mat<real>>(
			ctx.settings.defaultTolerance,
			mat_distance<mat<real>>
		);

		ctx.equals("jacobian(identity_2D)",
			J,
			expected,
			opt
		);
	}

	// Test Jacobian of scaled field [2x, 2y, 2z]: J = 2*I
	{
		real scalar = 2;

		auto f = [scalar](dvec v) -> dvec {
			return v * scalar;
		};

		vec<real> point = {
			1, 1, 1
		};

		auto J = jacobian(f, point);

		mat<real> expected = algebra::identity<mat<real>>(3, 3) * scalar;

		auto opt = prec::equation_options<mat<real>>(
			ctx.settings.defaultTolerance,
			mat_distance<mat<real>>
		);

		ctx.equals(
			"jacobian(scaled_field_3D)",
			J, expected, opt
		);
	}

	// Test Jacobian of f = [x^2, y^2]
	{
		auto f = [](dvec v) -> dvec {
			return { v[0] * v[0], v[1] * v[1] };
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto J = jacobian(f, point);

		mat<real> expected(2, 2);
		expected(0, 0) = 2.0 * point[0]; expected(0, 1) = 0.0;
		expected(1, 0) = 0.0; expected(1, 1) = 2.0 * point[1];

		auto opt = prec::equation_options<mat<real>>(
			ctx.settings.defaultTolerance,
			mat_distance<mat<real>>
		);

		ctx.equals("jacobian([x2,y2], [2,3])",
			J,
			expected,
			opt
		);
	}

	// Test Jacobian with constant component: f = [0, y, -x]
	{
		auto f = [](dvec v) -> dvec {
			return { 0.0, v[1], -v[0] };
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto J = jacobian(f, point);

		mat<real> expected(3, 3);
		expected(0, 0) = 0.0; expected(0, 1) = 0.0; expected(0, 2) = 0.0;
		expected(1, 0) = 0.0; expected(1, 1) = 1.0; expected(1, 2) = 0.0;
		expected(2, 0) = -1.0; expected(2, 1) = 0.0; expected(2, 2) = 0.0;

		auto opt = prec::equation_options<mat<real>>(
			ctx.settings.defaultTolerance,
			mat_distance<mat<real>>
		);

		ctx.equals("jacobian([0,y,-x])",
			J,
			expected,
			opt
		);
	}

	// Test curl of conservative field: curl(f) = 0
	{
		auto f = [](dvec v) -> dvec {
			return 2 * v;
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto c = curl(f, point);
		vec<real> expected = {0.0, 0.0, 0.0};

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);

		ctx.equals(
			"curl(f) = 0",
			c, expected, opt
		);
	}

	// Test curl of [x, y, z]: curl([x,y,z]) = 0
	{
		auto f = [](dvec v) { return v; };

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};
		
		auto c = curl(f, point);

		vec<real> expected = {0.0, 0.0, 0.0};

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);

		ctx.equals("curl([x,y,z])", c, expected, opt);
	}

	// Test curl of [-y, x, 0]: curl([-y,x,0]) = [0, 0, 2]
	{
		auto f = [](dvec v) -> dvec {

			return {
				-v[1], v[0], dreal(0.0)
			};
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto c = curl(f, point);
		vec<real> expected = {0.0, 0.0, 2.0};

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);
		
		ctx.equals("curl([-y,x,0])",c, expected, opt);
	}

	// Test curl of [0, z, -y]: curl([0,z,-y]) = [-2, 0, 0]
	{
		auto f = [](dvec v) -> dvec { 
			
			return {
				0.0, v[2], -v[1]
			};
		};

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		auto c = curl(f, point);
		vec<real> expected = {-2.0, 0.0, 0.0};

		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);
		
		ctx.equals("curl([0,z,-y])", c, expected, opt);
	}

	
	// dual.h

	// Test basic dual arithmetic
	{
		real a_re = rnd.gaussian(0, 1E+06);
		real a_dual = rnd.gaussian(0, 1E+06);
		real b_re = rnd.gaussian(0, 1E+06);
		real b_dual = rnd.gaussian(0, 1E+06);

		dual a(a_re, a_dual);
		dual b(b_re, b_dual);
		
		ctx.equals("dual::operator+()", a.Re() + b.Re(), a_re + b_re); 
		ctx.equals("dual::operator*()", (a * b).Re(), a_re * b_re);
		ctx.equals("dual::operator/()", (a / b).Re(), a_re / b_re);
	}

	// Test deriv() with polynomial
	{
		auto f = [](dual x) { 
			return x * x - 3 * x + 2.0; 
		};

		real x1 = rnd.gaussian(0, MAX);
		ctx.equals("deriv(poly2) (1)", 
			deriv(f, x1),
			dpoly2(x1)
		);

		auto f3 = [](dual x) {
			return x * x * x + 2 * x * x - x + dual(5.0);
		};

		real x2 = rnd.gaussian(0, MAX);
		ctx.equals("deriv(poly3)",
			deriv(f3, x2),
			dpoly3(x2)
		);
	}

	// Test deriv() with trigonometric functions
	{
		auto sin_x = [](dual x) { return sin(x); };
		auto cos_x = [](dual x) { return cos(x); };

		real x_sin = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv(sin, random)",
			deriv(sin_x, x_sin),
			std::cos(x_sin)
		);

		real x_cos = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv(cos, random)",
			deriv(cos_x, x_cos),
			-std::sin(x_cos)
		);

		real x_sin2 = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv(sin, random2)",
			deriv(sin_x, x_sin2),
			std::cos(x_sin2)
		);
	}

	// Test deriv() with exponential and logarithm
	{
		auto exp_x = [](dual x) { return exp(x); };
		auto ln_x = [](dual x) { return ln(x); };

		real x_exp1 = rnd.uniform(-10.0, 10.0);
		ctx.equals("deriv(exp, random)",
			deriv(exp_x, x_exp1),
			std::exp(x_exp1)
		);

		real x_exp2 = rnd.uniform(0.1, 4.1);
		ctx.equals("deriv(exp, random2)",
			deriv(exp_x, x_exp2),
			std::exp(x_exp2)
		);

		real x_ln = rnd.uniform(0.5, 10.0);
		ctx.equals("deriv(ln, random)",
			deriv(ln_x, x_ln),
			1.0 / x_ln
		);
	}

	// Test deriv() with sqrt
	{
		auto sqrt_x = [](dual x) { return sqrt(x); };

		real x_sqrt1 = rnd.uniform(0.1, 1E+06);
		ctx.equals("deriv(sqrt, random)",
			deriv(sqrt_x, x_sqrt1),
			0.5 / std::sqrt(x_sqrt1)
		);

		real x_sqrt2 = rnd.uniform(0.1, 1E+06);
		ctx.equals("deriv(sqrt, random2)",
			deriv(sqrt_x, x_sqrt2),
			0.5 / std::sqrt(x_sqrt2)
		);
	}

	// Test deriv() with product rule
	{
		auto f = [](dual x) { return x * sin(x); };

		real x = rnd.uniform(-1E+06, +1E+06);
		real expected = std::sin(x) + x * std::cos(x);

		ctx.equals("deriv(x * sin(x), random)",
			deriv(f, x),
			expected
		);
	}

	// Test deriv() with quotient rule
	{
		auto f = [](dual x) { return sin(x) / (1 + cos(x)); };

		real x = rnd.uniform(0.1, 5.9);
		real sin_x = std::sin(x);
		real cos_x = std::cos(x);
		real expected = (cos_x * (1 + cos_x) - sin_x * (-sin_x)) / ((1 + cos_x) * (1 + cos_x));

		ctx.equals("deriv(sin(x)/(1+cos(x)), random)",
			deriv(f, x),
			expected
		);
	}

	// Test deriv() with chain rule
	{
		auto f = [](dual x) {
			return sin(x * x);
		};

		real x = rnd.uniform(-2.0, 2.0);
		real expected = 2 * x * std::cos(x * x);

		ctx.equals("deriv(sin(x^2), random)",
			deriv(f, x),
			expected
		);
	}

	{
		dual d(3.0, 4.0);
		dual e(-2.0, 5.0);

		ctx.equals("dual::operator==", d == dual(3.0, 4.0), true);
		ctx.equals("dual::operator!=", d == e, false);
		ctx.equals("dual::Re()", d.Re(), 3.0);
		ctx.equals("dual::Dual()", d.Dual(), 4.0);

		dual access;
		access.Re() = rnd.uniform(-5.0, 5.0);
		access.Dual() = rnd.uniform(-5.0, 5.0);
		ctx.equals("dual::Re()", access.Re(), access.a);
		ctx.equals("dual::Dual()", access.Dual(), access.b);

		ctx.equals("conjugate(dual)", conjugate(d).Dual(), -4.0);
		ctx.equals("dual::inverse() real", d.inverse().Re(), 1.0 / 3.0);
		ctx.equals("dual::inverse() dual", d.inverse().Dual(), -4.0 / 9.0);
		ctx.equals("dual::operator-", (-d).Re(), -3.0);
		ctx.equals("dual::operator+", (d + e).Re(), 1.0);
		ctx.equals("dual::operator-", (d - e).Dual(), -1.0);
		ctx.equals("dual::operator*", (d * 2.0).Dual(), 8.0);
		ctx.equals("dual::operator/", (d / 2.0).Dual(), 2.0);
		ctx.equals("real + dual", (2.0 + d).Re(), 5.0);
		ctx.equals("real - dual", (2.0 - d).Re(), -1.0);

		vec2 v = d.to_vec();
		dual from_v;
		from_v.from_vec(v);
		ctx.equals("dual::to_vec()/from_vec() real", from_v.Re(), 3.0);
		ctx.equals("dual::to_vec()/from_vec() dual", from_v.Dual(), 4.0);

		mat2 m = d.to_mat();
		ctx.equals("dual::to_mat() 00", m(0, 0), 3.0);
		ctx.equals("dual::to_mat() 01", m(0, 1), 4.0);
		ctx.equals("dual::to_mat() 10", m(1, 0), 0.0);
		ctx.equals("dual::to_mat() 11", m(1, 1), 3.0);
	}

	// dual_functions.h coverage
	{
		real x = rnd.uniform(0.5, 3.0);
		dual dx(x, 1.0);

		ctx.equals("square(dual)", square(dx).Re(), x * x);
		ctx.equals("cube(dual)", cube(dx).Re(), x * x * x);
		ctx.equals("pow(dual,5)", pow(dx, 5).Re(), th::pow(x, 5.0));
		ctx.equals("pow(dual,5) d", pow(dx, 5).Dual(), 5.0 * th::pow(x, 4.0));
		ctx.equals("log2(dual)", log2(dx).Dual(), LOG2E / x);
		ctx.equals("log10(dual)", log10(dx).Dual(), LOG10E / x);
		ctx.equals("sinh(dual)", sinh(dx).Dual(), th::cosh(x));
		ctx.equals("cosh(dual)", cosh(dx).Dual(), th::sinh(x));
		ctx.equals("tanh(dual)", tanh(dx).Dual(), 1.0 / square(th::cosh(x)));
	}

	{
		real x = rnd.uniform(0.2, 0.8);
		dual dx(x, 1.0);

		ctx.equals("tan(dual)", tan(dx).Dual(), 1.0 / square(th::cos(x)));
		ctx.equals("cot(dual)", cot(dx).Dual(), -1.0 / square(th::sin(x)));
		ctx.equals("asin(dual)", asin(dx).Dual(), 1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("acos(dual)", acos(dx).Dual(), -1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("atan(dual)", atan(dx).Dual(), 1.0 / (1.0 + square(x)));
	}

	{
		real x = rnd.gaussian(0, MAX);

		if (th::abs(x) < MACH_EPSILON)
			x = 1.0;
		
		dual dx(x, 1.0);
		ctx.equals("abs(dual)", abs(dx).Dual(), x > 0 ? 1.0 : -1.0);
	}

	
	// dual2.h

	// Test deriv2() with polynomial
	{
		auto f = [](dual2 x) { 
			return x * x * x + 2 * x * x - x + 5.0;
		};

		real x1 = rnd.gaussian(0, MAX);
		ctx.equals("deriv2(poly3, random)",
			deriv2(f, x1),
			d2poly3(x1)
		);

		real x2 = rnd.gaussian(0, MAX);
		ctx.equals("deriv2(poly3, random2)",
			deriv2(f, x2),
			d2poly3(x2)
		);
	}

	// Test deriv2() with trigonometric functions
	{
		auto sin_x = [](dual2 x) { return sin(x); };
		auto cos_x = [](dual2 x) { return cos(x); };

		real x_sin = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv2(sin, random)",
			deriv2(sin_x, x_sin),
			-std::sin(x_sin)
		);

		real x_sin2 = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv2(sin, random2)",
			deriv2(sin_x, x_sin2),
			-std::sin(x_sin2)
		);

		real x_cos = rnd.uniform(-th::PI, th::PI);
		ctx.equals("deriv2(cos, random)",
			deriv2(cos_x, x_cos),
			-std::cos(x_cos)
		);
	}

	// Test deriv2() with exponential
	{
		auto exp_x = [](dual2 x) { return exp(x); };

		real x = rnd.uniform(-100.0, 10.0);

		ctx.equals("deriv2(exp, random)",
			deriv2(exp_x, x),
			std::exp(x)
		);
	}

	// Test deriv2() with logarithm
	{
		auto ln_x = [](dual2 x) { return ln(x); };

		real x = rnd.uniform(0.5, 10.0);

		ctx.equals("deriv2(ln, random)",
			deriv2(ln_x, x),
			-1.0 / (x * x)
		);
	}

	// Test deriv2() with product x * sin(x)
	{
		auto f = [](dual2 x) { return x * sin(x); };

		real x = rnd.uniform(-th::PI, th::PI);
		real expected = 2 * std::cos(x) - x * std::sin(x);

		ctx.equals("deriv2(x*sin(x), random)",
			deriv2(f, x),
			expected
		);
	}

	{
		real x = rnd.uniform(0.5, MAX);
		dual2 dx(x, 1.0, 0.0);

		ctx.equals("dual2::operator==()", dx == dual2(x, 1.0, 0.0), true);
		ctx.equals("dual2::Re()", dx.Re(), x);
		ctx.equals("dual2::Dual1()", dx.Dual1(), 1.0);
		ctx.equals("dual2::Dual2()", dx.Dual2(), 0.0);
		ctx.equals("dual2::conjugate().Dual1()", dx.conjugate().Dual1(), -1.0);
		ctx.equals("dual2::inverse().Re()", dx.inverse().Re(), 1.0 / x);

		vec3 v = dx.to_vec();
		dual2 from_v;
		from_v.from_vec(v);

		ctx.equals("dual2::to_vec() / ::from_vec() real", from_v.Re(), x);
		ctx.equals("dual2::to_vec() / ::from_vec() dual1", from_v.Dual1(), 1.0);
		ctx.equals("dual2::to_vec() / ::from_vec() dual2", from_v.Dual2(), 0.0);
	}

	{
		real x = rnd.uniform(0.5, MAX);
		dual2 dx(x, 1.0, 0.0);

		ctx.equals("square(dual2) d2", square(dx).Dual2(), 2.0);
		ctx.equals("cube(dual2) d2", cube(dx).Dual2(), 6.0 * x);
		ctx.equals("pow(dual2,5) d2", pow(dx, 5).Dual2(), 20.0 * th::pow(x, 3.0));
		ctx.equals("ln(dual2) d1", ln(dx).Dual1(), 1.0 / x);
		ctx.equals("ln(dual2) d2", ln(dx).Dual2(), -1.0 / square(x));
		ctx.equals("log2(dual2) d1", log2(dx).Dual1(), LOG2E / x);
		ctx.equals("log10(dual2) d1", log10(dx).Dual1(), LOG10E / x);

		x = rnd.uniform(0.5, 10.0);
		dx = dual2(x, 1.0, 0.0);

		ctx.equals("exp(dual2) d1", exp(dx).Dual1(), th::exp(x));
		ctx.equals("exp(dual2) d2", exp(dx).Dual2(), th::exp(x));
		ctx.equals("sinh(dual2) d1", sinh(dx).Dual1(), th::cosh(x));
		ctx.equals("cosh(dual2) d1", cosh(dx).Dual1(), th::sinh(x));
		ctx.equals("tanh(dual2) d1", tanh(dx).Dual1(), 1.0 / square(th::cosh(x)));
	}

	{
		real x = rnd.uniform(0.1, 0.9);
		dual2 dx(x, 1.0, 0.0);

		ctx.equals("tan(dual2) d1", tan(dx).Dual1(), 1.0 / square(th::cos(x)));
		ctx.equals("cot(dual2) d1", cot(dx).Dual1(), -1.0 / square(th::sin(x)));
		ctx.equals("asin(dual2) d1", asin(dx).Dual1(), 1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("acos(dual2) d1", acos(dx).Dual1(), -1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("atan(dual2) d1", atan(dx).Dual1(), 1.0 / (1.0 + square(x)));
	}

	{
		real x = rnd.gaussian(0, MAX);

		if (th::abs(x) < MACH_EPSILON)
			x = 1.0;
		
		dual2 dx(x, 1.0, 0.0);
		ctx.equals("abs(dual2) d1", abs(dx).Dual1(), x > 0 ? 1.0 : -1.0);
		ctx.equals("abs(dual2) d2", abs(dx).Dual2(), 0.0);
	}

	// multidual.h and multidual_functions.h

	{
		real x = rnd.uniform(0.5, 10.0);

		vec<real, 2> basis = {1.0, 0.0};
		multidual<2> dx(x, basis);

		ctx.equals("square(multidual)", square(dx).Dual(0), 2.0 * x);
		ctx.equals("cube(multidual)", cube(dx).Dual(0), 3.0 * x * x);
		ctx.equals("pow(multidual,5)", pow(dx, 5).Dual(0), 5.0 * th::pow(x, 4.0));
		ctx.equals("exp(multidual)", exp(dx).Dual(0), th::exp(x));
		ctx.equals("ln(multidual)", ln(dx).Dual(0), 1.0 / x);
		ctx.equals("log2(multidual)", log2(dx).Dual(0), LOG2E / x);
		ctx.equals("log10(multidual)", log10(dx).Dual(0), LOG10E / x);
		ctx.equals("sinh(multidual)", sinh(dx).Dual(0), th::cosh(x));
		ctx.equals("cosh(multidual)", cosh(dx).Dual(0), th::sinh(x));
		ctx.equals("tanh(multidual)", tanh(dx).Dual(0), 1.0 / square(th::cosh(x)));
	}

	{
		real x = rnd.uniform(0.1, 0.9);

		vec<real, 2> basis = {1.0, 0.0};
		multidual<2> dx(x, basis);

		ctx.equals("tan(multidual)", tan(dx).Dual(0), 1.0 / square(th::cos(x)));
		ctx.equals("cot(multidual)", cot(dx).Dual(0), -1.0 / square(th::sin(x)));
		ctx.equals("asin(multidual)", asin(dx).Dual(0), 1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("acos(multidual)", acos(dx).Dual(0), -1.0 / th::sqrt(1.0 - square(x)));
		ctx.equals("atan(multidual)", atan(dx).Dual(0), 1.0 / (1.0 + square(x)));
	}

	{
		real x = rnd.gaussian(0, MAX);

		if (th::abs(x) < MACH_EPSILON)
			x = 1.0;
		
		vec<real, 2> basis = {1.0, 0.0};
		multidual<2> dx(x, basis);

		ctx.equals("abs(multidual)", abs(dx).Dual(0), x > 0 ? 1.0 : -1.0);
	}

	// Regression test: scalar multiplication must preserve dual vector values
	{
		vec<real, 3> basis = {1.0, -2.0, 3.0};
		multidual<3> dx(4.0, basis);

		multidual<3> left = dx * 2.0;
		multidual<3> right = 2.0 * dx;

		ctx.equals("multidual * real real", left.Re(), 8.0);
		ctx.equals("multidual * real dual0", left.Dual(0), 2.0);
		ctx.equals("multidual * real dual1", left.Dual(1), -4.0);
		ctx.equals("multidual * real dual2", left.Dual(2), 6.0);
		ctx.equals("real * multidual real", right.Re(), 8.0);
		ctx.equals("real * multidual dual0", right.Dual(0), 2.0);
		ctx.equals("real * multidual dual1", right.Dual(1), -4.0);
		ctx.equals("real * multidual dual2", right.Dual(2), 6.0);
	}

	// Regression test: inverse and division should agree on both parts
	{
		vec<real, 2> basis = {2.0, -1.0};
		multidual<2> dx (5.0, basis);

		multidual<2> inv_dx = dx.inverse();
		multidual<2> div_dx = dx / 2.0;
		multidual<2> div_real = 10.0 / dx;

		ctx.equals("multidual::inverse().Re()", inv_dx.Re(), 0.2);
		ctx.equals("multidual::inverse().Dual(0)", inv_dx.Dual(0), -0.08);
		ctx.equals("multidual::inverse().Dual(1)", inv_dx.Dual(1), 0.04);

		ctx.equals("multidual::operator/(real).Re()", div_dx.Re(), 2.5);
		ctx.equals("multidual::operator/(real).Dual(0)", div_dx.Dual(0), 1.0);
		ctx.equals("multidual::operator/(real).Dual(1)", div_dx.Dual(1), -0.5);

		ctx.equals("multidual::operator/(multidual).Re()", div_real.Re(), 2.0);
		ctx.equals("multidual::operator/(multidual).Dual(0)", div_real.Dual(0), -0.8);
		ctx.equals("multidual::operator/(multidual).Dual(1)", div_real.Dual(1), 0.4);
	}

	// tanh(dual) should stay accurate for large inputs
	{
		const real x = rnd.gaussian(0, MAX);

		dual dx (x, 1);
		dual y = tanh(dx);

		const real expected_re = th::tanh(x);
		const real expected_dual = 1.0 / square(th::cosh(x));

		ctx.equals("tanh(dual) (stability, real)", y.Re(), expected_re);
		ctx.equals("tanh(dual) (stability, dual)", y.Dual(), expected_dual);
	}

	// tanh(multidual) should stay accurate for large inputs
	{
		const real x = rnd.gaussian(0, MAX);

		vec<real, 1> basis = {1.0};
		multidual<1> dx (x, basis);
		multidual<1> y = tanh(dx);

		const real expected_re = th::tanh(x);
		const real expected_dual = 1.0 / square(th::cosh(x));

		ctx.equals("tanh(multidual) (stability, real)", y.Re(), expected_re);
		ctx.equals("tanh(multidual) (stability, dual)", y.Dual(0), expected_dual);
	}

	// Currying overloads
	{
		auto f = [](dual x) { return x * x * x + dual(2.0) * x; };

		auto df = deriv(f);
		real x = rnd.uniform(-MAX, MAX);
		
		ctx.equals("deriv(f)", df(x), 3.0 * x * x + 2.0);
	}

	{
		auto f = [](dual2 x) { return x * x * x + dual2(2.0) * x; };
		auto d2f = deriv2(f);
		
		real x = rnd.uniform(-MAX, MAX);

		ctx.equals("deriv2(f)", d2f(x), 6.0 * x, 1e-10);
	}

	{
		auto f = [](dvec v) { return v[0] * v[0] + v[1] * v[1]; };
		auto grad_f = gradient(f);

		vec<real> point = { rnd.uniform(-4.0, 4.0), rnd.uniform(-4.0, 4.0) };
		vec<real> expected = { 2.0 * point[0], 2.0 * point[1] };
		
		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);
		
		ctx.equals("gradient(f)", grad_f(point), expected, opt);
	}

	{
		auto f = [](dvec v) -> dvec {
			return v;
		};

		auto Jf = jacobian(f);
		
		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		mat<real> expected(3, 3);
		expected(0, 0) = 1.0; expected(0, 1) = 0.0; expected(0, 2) = 0.0;
		expected(1, 0) = 0.0; expected(1, 1) = 1.0; expected(1, 2) = 0.0;
		expected(2, 0) = 0.0; expected(2, 1) = 0.0; expected(2, 2) = 1.0;

		auto opt = prec::equation_options<mat<real>>(
			ctx.settings.defaultTolerance,
			mat_distance<mat<real>>
		);

		ctx.equals("jacobian(f)", Jf(point), expected, opt);
	}

	{
		auto f = [](dvec v) -> dvec {

			return {
				-v[1], v[0], dreal(0.0)
			};
		};

		auto cf = curl(f);
		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};
		
		vec<real> expected = {0.0, 0.0, 2.0};
		auto opt = prec::equation_options<vec<real>>(
			ctx.settings.defaultTolerance,
			prec::distance::euclidean<vec<real>>
		);

		ctx.equals("curl(f)", cf(point), expected, opt);
	}

	{
		auto f = [](vec<dual2> v) {
			return v[0] * v[0] + v[1] * v[1];
		};
		
		auto lap_f = laplacian(f);

		vec<real> point = {
			rnd.gaussian(0, MAX),
			rnd.gaussian(0, MAX)
		};

		ctx.equals("laplacian(f)", lap_f(point), 2.0 * point.size());
	}


	// Precision estimates for autodiff operators
	{
		auto f = [](dual x) {
			return x * x * x + 2 * x * x - x + dual(5.0);
		};

		auto opt = prec::estimate_options<real, real>(
			prec::interval(-MAX, MAX),
			prec::estimator::quadrature1D()
		);

		ctx.estimate("deriv(poly3)", deriv(f), dpoly3, opt);
	}

	{
		auto f = [](dual2 x) {
			return x * x * x + 2 * x * x - x + dual2(5.0);
		};

		auto opt = prec::estimate_options<real, real>(
			prec::interval(-10.0, 10.0),
			prec::estimator::quadrature1D()
		);

		ctx.estimate("deriv2(poly3)", deriv2(f), d2poly3, opt);
	}

	{
		auto f = [](dvec v) {
			return v[0] * v[0] + v[1] * v[1];
		};

		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(-MAX, MAX), prec::interval(-MAX, MAX) },
			prec::estimator::montecarlo<real>(ctx.random, 2)
		);

		ctx.homogeneous(
			"gradient(x2+y2) residual",
			[f](std::vector<real> x) {
				return algebra::norm(
					autodiff::gradient(f, vec<real>(x)) - vec<real>(2.0 * x[0], 2.0 * x[1])
				);
			},
			opt
		);
	}

	{
		auto f = [](dvec v) -> dvec {
			return { -v[1], v[0], dreal(0.0) };
		};

		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(-MAX, MAX), prec::interval(-MAX, MAX), prec::interval(-MAX, MAX) },
			prec::estimator::montecarlo<real>(ctx.random, 3)
		);

		ctx.homogeneous(
			"divergence([-y,x,0]) residual",
			[f](std::vector<real> x) {
				return divergence(f, vec<real>(x));
			},
			opt
		);
	}

	{
		auto f = [](dvec v) -> dvec {
			return { 2.0 * v[0], 2.0 * v[1], 2.0 * v[2] };
		};

		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(-MAX, MAX), prec::interval(-MAX, MAX), prec::interval(-MAX, MAX) },
			prec::estimator::montecarlo<real>(ctx.random, 3)
		);

		ctx.homogeneous(
			"curl(grad field) residual",
			[f](std::vector<real> x) {
				return curl(f, vec<real>(x)).norm();
			},
			opt
		);
	}

	{
		auto f = [](vec<dual2> v) -> dual2 {
			return ln(v[0] * v[0] + v[1] * v[1]);
		};

		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(0.5, MAX), prec::interval(0.5, MAX) },
			prec::estimator::montecarlo<real>(ctx.random, 2)
		);

		ctx.homogeneous(
			"laplacian(ln(r2)) residual",
			[f](std::vector<real> x) {
				return laplacian(f, vec<real>(x));
			},
			opt
		);
	}
}
