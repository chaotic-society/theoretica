
///
/// @file test_polynomial.cpp Polynomial class test cases
///

#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


polynomial<real> rand_polyn(random::random_source& rnd, unsigned int degree, real stdev = 1E+06) {

	std::vector<real> coeffs (degree + 1);
	for (unsigned int i = 0; i <= degree; ++i)
		coeffs[i] = rnd.gaussian(0, stdev);

	return polynomial<real>(std::move(coeffs));
}


// Distance function for polynomial comparison
long double distance_polyn(const polynomial<real>& p1, const polynomial<real>& p2) {

	const polynomial<real> d = p1 - p2;
	real r = -inf();

	for (size_t i = 0; i < d.size(); ++i)
		r = max(r, th::abs(d[i]));

	return r;
}


auto estimate_polyn(random::random_context& rctx, unsigned int degree, real stdev = 1E+06) {

	return [&rctx, degree, stdev](
		std::function<real(polynomial<real>)> f, std::function<real(polynomial<real>)> f_exact,
		prec::estimate_options<real, polynomial<real>> opt) {

		auto rnd = rctx.get_rnd();
		prec::estimate_result res;

		real sum = 0;
		real sum_sqr = 0;
		real max = 0;

		for (unsigned int i = 0; i < opt.iterations; ++i) {

			polynomial<real> p = rand_polyn(rnd, degree, stdev);

			real val = std::abs(f(p) - f_exact(p));
			sum += val;
			sum_sqr += val * val;
			max = std::max(max, val);
		}

		res.meanErr = sum / opt.iterations;
		res.rmsErr = std::sqrt(sum_sqr / opt.iterations);
		res.maxErr = max;
		res.failed = opt.fail(res);

		return res;
	};
}


auto make_options(random::random_context& rctx, unsigned int degree, real stdev = 1E+06) {

	return prec::estimate_options<real, polynomial<real>>(
		prec::interval(-stdev, +stdev),
		estimate_polyn(rctx, degree, stdev),
		1E-08,
		1000
	);
}


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("polynomial", argc, argv);
	ctx.settings.defaultIterations = 1000;
	ctx.output->settings.outputFiles = { "test/prec/prec_polynomial.csv" };
	
	random::random_source rnd = ctx.random->get_rnd();
	const real MAX = 1E+06;
	auto zero = [](polynomial<real> p) -> real { return 0.0; };

	auto polyn_opt = prec::equation_options<polynomial<real>>(
		1E-08, distance_polyn
	);

	// polynomial.h

	// Default constructor
	{
		polynomial<real> p;
		ctx.equals("polynomial()", p.size(), 1);
		ctx.equals("polynomial()  value", p[0], 0.0);
	}

	// Constant constructor
	{
		real c = rnd.gaussian(0, MAX);
		polynomial<real> p (c);
		ctx.equals("polynomial() (size)", p.size(), 1);
		ctx.equals("polynomial() (value)", p[0], c);
	}

	// Vector constructor
	{
		vec<real> coeffs = {rnd.gaussian(0, MAX), rnd.gaussian(0, MAX), rnd.gaussian(0, MAX)};
		polynomial<real> p(coeffs);

		ctx.equals("polynomial(coeffs) (size)", p.size(), 3);
		ctx.equals("polynomial(coeffs) [0]", p[0], coeffs[0]);
		ctx.equals("polynomial(coeffs) [1]", p[1], coeffs[1]);
		ctx.equals("polynomial(coeffs) [2]", p[2], coeffs[2]);
	}

	// Initializer list constructor
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real c3 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2, c3};

		ctx.equals("polynomial(initializer_list) (size)", p.size(), 4);
		ctx.equals("polynomial(initializer_list) [0]", p[0], c0);
		ctx.equals("polynomial(initializer_list) [3]", p[3], c3);
	}


	// Element Access

	// operator[] read access
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};

		ctx.equals("polynomial::operator[0]", p[0], c0);
		ctx.equals("polynomial::operator[1]", p[1], c1);
		ctx.equals("polynomial::operator[2]", p[2], c2);
	}

	// operator[] write access
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real new_val = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		p[1] = new_val;
		ctx.equals("polynomial::operator[]", p[1], new_val);
	}

	// at() read access
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		ctx.equals("polynomial::at(0)", p.at(0), c0);
		ctx.equals("polynomial::at(2)", p.at(2), c2);
	}

	// at() write access
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real new_val = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		p.at(1) = new_val;
		ctx.equals("polynomial::at()", p.at(1), new_val);
	}


	// Evaluation

	// eval() constant polynomial
	{
		real c = rnd.gaussian(0, MAX);
		real x = rnd.gaussian(0, MAX);

		polynomial<real> p = {c};
		ctx.equals("eval(const)", p.eval(x), c);
	}

	// eval() at zero
	{
		ctx.estimate("eval(0)", [](polynomial<real> p) {
			return real(p.eval(0.0) - p[0]);
		}, zero, make_options(*ctx.random, 10));
	}

	// eval() linear polynomial: c0 + c1*x
	{
		real c0 = rnd.gaussian(0, MAX), c1 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1};

		ctx.equals("eval(0) (linear)", p.eval(0.0), c0);
		ctx.equals("eval(1) (linear)", p.eval(1.0), c0 + c1);
		ctx.equals("eval(2) (linear)", p.eval(2.0), c0 + 2*c1);
		ctx.equals("eval(-1) (linear)", p.eval(-1.0), c0 - c1);
	}

	// eval() quadratic: c0 + c1*x + c2*x^2
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		
		polynomial<real> p = {c0, c1, c2};
		ctx.equals("eval(0) (quadratic)", p.eval(0.0), c0);
		ctx.equals("eval(1) (quadratic)", p.eval(1.0), c0 + c1 + c2);
		ctx.equals("eval(2) (quadratic)", p.eval(2.0), c0 + 2*c1 + 4*c2);
		ctx.equals("eval(-1) (quadratic)", p.eval(-1.0), c0 - c1 + c2);
	}

	// eval() cubic: c0 + c1*x + c2*x^2 + c3*x^3
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real c3 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2, c3};
		ctx.equals("eval(0) (cubic)", p.eval(0.0), c0);
		ctx.equals("eval(1) (cubic)", p.eval(1.0), c0 + c1 + c2 + c3);
		ctx.equals("eval(2) (cubic)", p.eval(2.0), c0 + 2*c1 + 4*c2 + 8*c3);
	}

	// operator() functional syntax
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		ctx.equals("operator(0.0)", p(0.0), c0);
		ctx.equals("operator(1.0)", p(1.0), c0 + c1 + c2);
		ctx.equals("operator(2.0)", p(2.0), c0 + 2*c1 + 4*c2);
	}

	// Properties

	// degree() constant
	{
		real c = rnd.gaussian(0, MAX);
		polynomial<real> p = {c};

		ctx.equals("degree(constant)", p.degree(), 0);
	}

	// degree() linear
	{
		real c0 = rnd.gaussian(0, MAX), c1 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1};
		ctx.equals("degree(linear)", p.degree(), 1);
	}

	// degree() quadratic
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};
		ctx.equals("degree(quadratic)", p.degree(), 2);
	}

	// degree() with trailing zeros
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2, 0.0, 0.0};

		ctx.equals("degree()", p.degree(), 2);
	}

	// trim() removes trailing zeros
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2, 0.0, 0.0};
		p.trim();
		ctx.equals("trim() (size)", p.size(), 3);
		ctx.equals("trim() (degree)", p.degree(), 2);
	}

	// trim() with all zeros except first
	{
		real c = rnd.gaussian(0, MAX);
		polynomial<real> p = {c, 0.0, 0.0};
		p.trim();
		ctx.equals("trim() (size)", p.size(), 1);
		ctx.equals("trim() (constant value)", p[0], c);
	}

	// size()
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real c3 = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2, c3};
		ctx.equals("size()", p.size(), 4);
	}

	// resize()
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		p.resize(5);
		ctx.equals("resize() (expand)", p.size(), 5);
		ctx.equals("resize() (expand, original)", p[0], c0);
	}

	// coeffs() const and non-const access
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real new_val = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		const polynomial<real>& p_const = p;

		ctx.equals("coeffs() const", p_const.coeffs()[1], c1);
		p.coeffs()[1] = new_val;
		ctx.equals("coeffs() non-const", p[1], new_val);
	}

	// Arithmetic Operations

	// Addition of two polynomials
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX), p1_2 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {p1_0, p1_1, p1_2};
		polynomial<real> p2 = {p2_0, p2_1};
		polynomial<real> result = p1 + p2;
		polynomial<real> expected = {p1_0 + p2_0, p1_1 + p2_1, p1_2};

		ctx.equals("operator+(polynomials)", result, expected, polyn_opt);
	}

	// Addition of different sizes
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX), p2_2 = rnd.gaussian(0, MAX), p2_3 = rnd.gaussian(0, MAX);
		polynomial<real> p1 = {p1_0, p1_1};
		polynomial<real> p2 = {p2_0, p2_1, p2_2, p2_3};
		polynomial<real> result = p1 + p2;
		ctx.equals("operator+() (different size, [0])", result[0], p1_0 + p2_0);
		ctx.equals("operator+() (different size, [3])", result[3], p2_3);
	}

	// Subtraction of two polynomials
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX), p1_2 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX);
		polynomial<real> p1 = {p1_0, p1_1, p1_2};
		polynomial<real> p2 = {p2_0, p2_1};
		polynomial<real> result = p1 - p2;
		polynomial<real> expected = {p1_0 - p2_0, p1_1 - p2_1, p1_2};
		ctx.equals("operator-(polynomials)", result, expected, polyn_opt);
	}

	// Multiplication of two polynomials
	{
		polynomial<real> p1 = {1.0, 1.0};  // 1 + x
		polynomial<real> p2 = {1.0, 1.0};  // 1 + x
		polynomial<real> result = p1 * p2;  // (1 + x)^2 = 1 + 2x + x^2
		polynomial<real> expected = {1.0, 2.0, 1.0};
		ctx.equals("operator*(polynomials)", result, expected, polyn_opt);
	}

	// Multiplication (1 + x) * (2 + 3x)
	{
		polynomial<real> p1 = {1.0, 1.0};
		polynomial<real> p2 = {2.0, 3.0};
		polynomial<real> result = p1 * p2;  // 2 + 3x + 2x + 3x^2 = 2 + 5x + 3x^2
		polynomial<real> expected = {2.0, 5.0, 3.0};
		ctx.equals("operator*(linear)", result, expected, polyn_opt);
	}

	// Scalar multiplication (left)
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real scalar = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};
		polynomial<real> result = p * scalar;
		polynomial<real> expected = {c0 * scalar, c1 * scalar, c2 * scalar};
		ctx.equals("operator*(scalar right)", result, expected, polyn_opt);
	}

	// Scalar multiplication (right, friend operator)
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real scalar = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};
		polynomial<real> result = scalar * p;
		polynomial<real> expected = {c0 * scalar, c1 * scalar, c2 * scalar};
		ctx.equals("operator*(scalar left)", result, expected, polyn_opt);
	}

	// Scalar division
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real scalar = rnd.gaussian(1, 1000);  // Avoid division by very small numbers
		polynomial<real> p = {c0, c1, c2};
		polynomial<real> result = p / scalar;
		polynomial<real> expected = {c0 / scalar, c1 / scalar, c2 / scalar};
		ctx.equals("operator/(scalar)", result, expected, polyn_opt);
	}

	// Polynomial division: (1 + 2x + x^2) / (1 + x) = 1 + x
	{
		polynomial<real> dividend = {1.0, 2.0, 1.0};  // 1 + 2x + x^2 = (1 + x)^2
		polynomial<real> divisor = {1.0, 1.0};        // 1 + x
		polynomial<real> result = dividend / divisor;
		polynomial<real> expected = {1.0, 1.0};
		ctx.equals("operator/(polynomial)", result, expected, polyn_opt);
	}

	// In-place Arithmetic

	// operator+=
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX), p2_2 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {p1_0, p1_1};
		polynomial<real> p2 = {p2_0, p2_1, p2_2};
		p1 += p2;

		polynomial<real> expected = {p1_0 + p2_0, p1_1 + p2_1, p2_2};
		ctx.equals("operator+=()", p1, expected, polyn_opt);
	}

	// operator-=
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX), p1_2 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {p1_0, p1_1, p1_2};
		polynomial<real> p2 = {p2_0, p2_1};
		p1 -= p2;
		polynomial<real> expected = {p1_0 - p2_0, p1_1 - p2_1, p1_2};
		ctx.equals("operator-=()", p1, expected, polyn_opt);
	}

	// operator*= (polynomial)
	{
		polynomial<real> p1 = {1.0, 1.0};  // 1 + x
		polynomial<real> p2 = {1.0, 1.0};  // 1 + x
		p1 *= p2;                          // (1 + x)^2 = 1 + 2x + x^2
		polynomial<real> expected = {1.0, 2.0, 1.0};
		ctx.equals("operator*=(polynomial)", p1, expected, polyn_opt);
	}

	// operator*= (scalar)
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real scalar = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};
		p *= scalar;
		polynomial<real> expected = {c0 * scalar, c1 * scalar, c2 * scalar};
		ctx.equals("operator*=(scalar)", p, expected, polyn_opt);
	}

	// operator/= (polynomial)
	{
		polynomial<real> p1 = {1.0, 2.0, 1.0};  // (1 + x)^2
		polynomial<real> p2 = {1.0, 1.0};       // 1 + x

		p1 /= p2;
		polynomial<real> expected = {1.0, 1.0};
		ctx.equals("operator/=(polynomial)", p1, expected, polyn_opt);
	}

	// operator/= (scalar)
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real scalar = rnd.uniform(1, 1E+03);  // Avoid division by very small numbers

		polynomial<real> p = {c0, c1, c2};
		p /= scalar;

		polynomial<real> expected = {c0 / scalar, c1 / scalar, c2 / scalar};
		ctx.equals("operator/=(scalar)", p, expected, polyn_opt);
	}

	// Equality

	// operator== equal polynomials
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		
		polynomial<real> p1 = {c0, c1, c2};
		polynomial<real> p2 = {c0, c1, c2};
		ctx.equals("operator==(equal)", p1 == p2, true);
	}

	// operator== unequal polynomials
	{
		real c0 = rnd.gaussian(0, MAX), c1 = rnd.gaussian(0, MAX), c2 = rnd.gaussian(0, MAX), c2_diff = rnd.gaussian(0, MAX);

		if (th::abs(c2_diff - c2) < 1E-06)
			c2_diff += 100;  // Ensure different

		polynomial<real> p1 = {c0, c1, c2};
		polynomial<real> p2 = {c0, c1, c2_diff};
		ctx.equals("operator==(unequal)", p1 == p2, false);
	}

	// operator== different sizes but equivalent
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {c0, c1, c2};
		polynomial<real> p2 = {c0, c1, c2, 0.0, 0.0};

		ctx.equals("operator==(trailing zeros)", p1 == p2, true);
	}

	// Root Finding

	// quadratic_roots: ax^2 + bx + c = 0
	{
		ctx.estimate("quadratic_roots()", [](polynomial<real> p) {

			auto roots = p.quadratic_roots();
			return max(
				th::abs(p(roots[0])),
				th::abs(p(roots[1]))
			);

		}, zero, make_options(*ctx.random, 2, 1E+03));
	}

	// quadratic_roots: x^2 + c = 0 (complex roots)
	{
		ctx.estimate("quadratic_roots() (complex)", [](polynomial<real> p) {

			p[0] = std::abs(p[0]);
			p[1] = 0.0;
			p[2] = 1.0;
			auto roots = p.quadratic_roots();

			return max(
				th::abs(p(roots[0])),
				th::abs(p(roots[1]))
			);

		}, zero, make_options(*ctx.random, 2, 1E+03));
	}

	// from_roots: construct from roots
	{
		real r1 = rnd.gaussian(0, 1e+03);
		real r2 = rnd.gaussian(0, 1e+03);
		vec<real> roots = {r1, r2};
		polynomial<real> p = polynomial<real>::from_roots(roots);
		
		// (x - r1)(x - r2) should have roots at r1 and r2
		ctx.equals("from_roots() (eval r1)", th::abs(p(r1)), 0.0);
		ctx.equals("from_roots() (eval r2)", th::abs(p(r2)), 0.0);

		// At x=0: p(0) = r1*r2
		ctx.equals("from_roots() (eval 0)", p(0.0), r1*r2);
	}

	// from_roots: multiple roots
	{
		real r1 = rnd.gaussian(0, 1e+02);
		real r2 = rnd.gaussian(0, 1e+02);
		real r3 = rnd.gaussian(0, 1e+02);

		vec<real> roots = {r1, r2, r3};
		polynomial<real> p = polynomial<real>::from_roots(roots);
		
		// (x - r1)(x - r2)(x - r3) should have roots at r1, r2, r3
		ctx.equals("from_roots(r1)", th::abs(p(r1)), 0.0);
		ctx.equals("from_roots(r2)", th::abs(p(r2)), 0.0);
		ctx.equals("from_roots(r3)", th::abs(p(r3)), 0.0);
	}

	// Static Methods

	// monomial
	{
		real coeff = rnd.gaussian(0, MAX);
		polynomial<real> p = polynomial<real>::monomial(coeff, 3);

		ctx.equals("monomial(size)", p.size(), 4);
		ctx.equals("monomial(value)", p[3], coeff);
		ctx.equals("monomial(other)", p[0], 0.0);
		ctx.equals("monomial(other2)", p[2], 0.0);
	}

	// monomial degree 0
	{
		real coeff = rnd.gaussian(0, MAX);
		polynomial<real> p = polynomial<real>::monomial(coeff, 0);
		ctx.equals("monomial(degree 0)", p[0], coeff);
	}

	// Iterators

	// begin/end iteration
	{
		polynomial<real> p = {1.0, 2.0, 3.0};
		int count = 0;
		for (auto it = p.begin(); it != p.end(); ++it) {
			count++;
		}
		ctx.equals("begin/end iteration count", count, p.size());
	}

	// String Conversion

	// to_string basic
	{
		polynomial<real> p = {1.0, 2.0};
		std::string s = p.to_string();
		ctx.equals("to_string (non-empty)", s.empty(), false);
	}


	// orthogonal.h - Legendre Polynomials

	// legendre_polynomial degree 0: P0(x) = 1
	{
		polynomial<real> P0 = legendre_polynomial(0);
		ctx.equals("legendre(0, at 0.5)", P0(0.5), 1.0, 1E-10);
		ctx.equals("legendre(0, at 1.0)", P0(1.0), 1.0, 1E-10);
	}

	// legendre_polynomial degree 1: P1(x) = x
	{
		polynomial<real> P1 = legendre_polynomial(1);
		ctx.equals("legendre(1, at 0.0)", P1(0.0), 0.0, 1E-10);
		ctx.equals("legendre(1, at 1.0)", P1(1.0), 1.0, 1E-10);
		ctx.equals("legendre(1, at 0.5)", P1(0.5), 0.5, 1E-10);
	}

	// legendre_polynomial degree 2: P2(x) = (3x^2 - 1)/2
	{
		polynomial<real> P2 = legendre_polynomial(2);

		// At x=0: P2(0) = -1/2
		ctx.equals("legendre(2, at 0)", P2(0.0), -0.5, 1E-10);

		// At x=1: P2(1) = 1
		ctx.equals("legendre(2, at 1)", P2(1.0), 1.0, 1E-10);
	}

	// legendre_polynomial degree 3
	{
		polynomial<real> P3 = legendre_polynomial(3);

		ctx.equals("legendre(3, size)", P3.degree(), 3);
	}

	// Laguerre polynomials

	// laguerre_polynomial degree 0: L0(x) = 1
	{
		polynomial<real> L0 = laguerre_polynomial(0);

		ctx.equals("laguerre(0, at 0.5)", L0(0.5), 1.0, 1E-10);
	}

	// laguerre_polynomial degree 1: L1(x) = 1 - x
	{
		polynomial<real> L1 = laguerre_polynomial(1);

		ctx.equals("laguerre(1, at 0)", L1(0.0), 1.0, 1E-10);
		ctx.equals("laguerre(1, at 1)", L1(1.0), 0.0, 1E-10);
	}

	// laguerre_polynomial degree 2
	{
		polynomial<real> L2 = laguerre_polynomial(2);

		ctx.equals("laguerre(2, at 0)", L2(0.0), 1.0, 1E-10);
	}

	// Hermite polynomials

	// hermite_polynomial degree 0: H0(x) = 1
	{
		polynomial<real> H0 = hermite_polynomial(0);

		ctx.equals("hermite(0, at 0.5)", H0(0.5), 1.0, 1E-10);
	}

	// hermite_polynomial degree 1: H1(x) = 2x
	{
		polynomial<real> H1 = hermite_polynomial(1);

		ctx.equals("hermite(1, at 0)", H1(0.0), 0.0, 1E-10);
		ctx.equals("hermite(1, at 1)", H1(1.0), 2.0, 1E-10);
	}

	// hermite_polynomial degree 2
	{
		polynomial<real> H2 = hermite_polynomial(2);

		ctx.equals("hermite(2, at 0)", H2(0.0), -2.0, 1E-10);
	}

	// Chebyshev polynomials (first kind)

	// chebyshev1_polynomial degree 0: T0(x) = 1
	{
		polynomial<real> T0 = chebyshev1_polynomial(0);

		ctx.equals("chebyshev1(0, at 0.5)", T0(0.5), 1.0, 1E-10);
	}

	// chebyshev1_polynomial degree 1: T1(x) = x
	{
		polynomial<real> T1 = chebyshev1_polynomial(1);

		ctx.equals("chebyshev1(1, at 0)", T1(0.0), 0.0, 1E-10);
		ctx.equals("chebyshev1(1, at 1)", T1(1.0), 1.0, 1E-10);
	}

	// chebyshev1_polynomial degree 2: T2(x) = 2x^2 - 1
	{
		polynomial<real> T2 = chebyshev1_polynomial(2);

		ctx.equals("chebyshev1(2, at 0)", T2(0.0), -1.0, 1E-10);
		ctx.equals("chebyshev1(2, at 1)", T2(1.0), 1.0, 1E-10);
	}

	// chebyshev1_polynomial property: |T_n(x)| <= 1 for |x| <= 1
	{
		polynomial<real> T3 = chebyshev1_polynomial(3);

		for (real x = -1.0; x <= 1.0; x += 0.2) {
			real val = T3(x);
			ctx.equals("chebyshev1(3) (bounded)", (val >= -1.0 && val <= 1.0), true);
		}
	}

	// Chebyshev polynomials (second kind)

	// chebyshev2_polynomial degree 0: U0(x) = 1
	{
		polynomial<real> U0 = chebyshev2_polynomial(0);

		ctx.equals("chebyshev2(0)", U0(0.5), 1.0, 1E-10);
	}

	// chebyshev2_polynomial degree 1: U1(x) = 2x
	{
		polynomial<real> U1 = chebyshev2_polynomial(1);

		ctx.equals("chebyshev2(1) (eval 0)", U1(0.0), 0.0, 1E-10);
		ctx.equals("chebyshev2(1) (eval 1)", U1(1.0), 2.0, 1E-10);
	}

	// general_laguerre_polynomial with alpha
	{
		polynomial<real> L0 = general_laguerre_polynomial(0.5, 0);

		ctx.equals("general_laguerre(alpha=0.5, n=0, at 0.5)", L0(0.5), 1.0, 1E-10);
	}

	// orthogonal.h - Legendre Roots and Weights

	// legendre_roots returns correct number of roots
	{
		auto roots = legendre_roots(3);

		ctx.equals("legendre_roots(3) (size)", roots.size(), 3);
	}

	// legendre_roots are within [-1, 1]
	{
		auto roots = legendre_roots(5);

		for (auto r : roots) {
			ctx.equals("legendre_roots(5) (bounded)", (r >= -1.0 && r <= 1.0), true);
		}
	}

	// legendre_weights returns correct number of weights
	{
		auto roots = legendre_roots(3);
		auto weights = legendre_weights(roots);

		ctx.equals("legendre_weights(3) (size)", weights.size(), 3);
	}

	// legendre_weights sum to 2
	{
		auto roots = legendre_roots(5);
		auto weights = legendre_weights(roots);
		
		real sum = 0.0;
		for (auto w : weights)
			sum += w;

		ctx.equals("legendre_weights(5) (sum)", sum, 2.0, 1E-10);
	}

	// hermite_weights returns correct number of weights
	{
		auto roots = legendre_roots(4);
		auto roots_h = legendre_roots(4);
		auto weights = hermite_weights(roots_h);

		ctx.equals("hermite_weights(4) (size)", weights.size(), roots_h.size());
	}

	// laguerre_weights returns correct number of weights
	{
		auto roots = legendre_roots(3);
		auto weights = laguerre_weights(roots);

		ctx.equals("laguerre_weights(3) (size)", weights.size(), roots.size());
	}


	// Edge Cases

	// eval() at x=0
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2};
		ctx.equals("eval(x=0)", p(0.0), c0);
	}

	// eval() high degree polynomial
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real c3 = rnd.gaussian(0, MAX);
		real c4 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2, c3, c4};
		real expected = c0 + 2*c1 + 4*c2 + 8*c3 + 16*c4;

		ctx.equals("eval() (high degree)", p(2.0), expected);
	}

	// Subtraction resulting in constant
	{
		real c0 = rnd.gaussian(0, MAX), c1 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {c0, c1};
		polynomial<real> p2 = {c0 - 2.0, c1};
		polynomial<real> result = p1 - p2;

		ctx.equals("operator-", result(1.0), 2.0);
	}

	// Multiple derivatives
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real c3 = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1, c2, c3};
		polynomial<real> dp = deriv(p);
		polynomial<real> d2p = deriv(dp);

		polynomial<real> expected_dp = {c1, 2*c2, 3*c3};
		polynomial<real> expected_d2p = {2*c2, 6*c3};

		ctx.equals("deriv(deriv)", d2p, expected_d2p, polyn_opt);
	}

	// polynomial construction with complex coefficients
	{
		real c0_r = rnd.gaussian(0, MAX), c0_i = rnd.gaussian(0, MAX);
		real c1_r = rnd.gaussian(0, MAX), c1_i = rnd.gaussian(0, MAX);

		polynomial<complex<>> p = {complex<>(c0_r, c0_i), complex<>(c1_r, c1_i)};
		ctx.equals("complex polynomial (size)", p.size(), 2);
	}

	// Polynomial addition with negative result
	{
		real p1_0 = rnd.gaussian(0, MAX), p1_1 = rnd.gaussian(0, MAX);
		real p2_0 = rnd.gaussian(0, MAX), p2_1 = rnd.gaussian(0, MAX);

		polynomial<real> p1 = {p1_0, p1_1};
		polynomial<real> p2 = {p2_0, p2_1};
		polynomial<real> result = p1 - p2;

		ctx.equals("operator- (negative)", result[0], p1_0 - p2_0);
	}

	// Trimming zero polynomial
	{
		polynomial<real> p = {0.0};
		p.trim();
		ctx.equals("trim() (zero poly)", p.size(), 1);
		ctx.equals("trim() (zero poly, value)", p[0], 0.0);
	}

	// Evaluation consistency: eval() vs operator()
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);
		real x = rnd.gaussian(0, MAX);
		polynomial<real> p = {c0, c1, c2};
		ctx.equals("eval vs operator()", p.eval(x), p(x));
	}

	// Constant polynomial after default construction
	{
		polynomial<real> p = polynomial<real>();
		ctx.equals("default construction (size)", p.size(), 1);
	}

	// Friend operator: scalar - polynomial
	{
		real c0 = rnd.gaussian(0, MAX), c1 = rnd.gaussian(0, MAX);
		real scalar = rnd.gaussian(0, MAX);

		polynomial<real> p = {c0, c1};
		polynomial<real> result = scalar - p;
		polynomial<real> expected = {scalar - c0, -c1};

		ctx.equals("friend operator-(scalar, poly)", result, expected, polyn_opt);
	}

	// polynomial with negative coefficients
	{
		real c0 = rnd.gaussian(0, MAX);
		real c1 = rnd.gaussian(0, MAX);
		real c2 = rnd.gaussian(0, MAX);

		polynomial<real> p = {-c0, -c1, -c2};
		real expected = -c0 - c1 - c2;

		ctx.equals("eval() (negative coefficients)", p(1.0), expected);
	}

	// Degree of zero polynomial
	{
		polynomial<real> p = polynomial<real>::monomial(0.0, 5);

		p.trim();
		ctx.equals("degree() (zero polynomial)", p.degree(), 0);
	}
	
}
