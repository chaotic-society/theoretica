
///
/// @file test_complex.cpp Test cases for complex numbers, quaternions, and phasors
///

#include "theoretica.h"
#include "chebyshev.h"
#include <cmath>

using namespace chebyshev;
using namespace theoretica;


const real MAX = 1E+09;
const real VARIANCE = 1E+04;


// Generate a random complex number
complex<real> rand_complex(random::random_source rnd, real variance = VARIANCE) {
	return complex<>(rnd.gaussian(0, variance), rnd.gaussian(0, variance));
}

// Generate a random phasor
phasor<real> rand_phasor(random::random_source rnd, real variance = VARIANCE) {
	return phasor<>(std::abs(rnd.gaussian(0, variance)), rnd.uniform(0, 2 * PI));
}

// Generate a random quaternion
quat<real> rand_quat(random::random_source rnd, real variance = VARIANCE) {
	return quat<>(rnd.gaussian(0, variance), rnd.gaussian(0, variance), rnd.gaussian(0, variance), rnd.gaussian(0, variance));
}


// Distance between two complex numbers
real_t distance_complex(complex<> z, complex<> w) {

	long double re = z.a - w.a;
	long double im = z.b - w.b;
	
	return std::sqrt(re*re + im*im);
}

// Distance between two phasors
real_t distance_phasor(phasor<> z, phasor<> w) {

	complex<> z1 (z.modulus * std::cos(z.phase), z.modulus * std::sin(z.phase));
	complex<> w1 (w.modulus * std::cos(w.phase), w.modulus * std::sin(w.phase));

	return distance_complex(z1, w1);
}

// Distance between two quaternions
real_t distance_quat(quat<> z, quat<> w) {

	long double diff1 = z.a - w.a;
	long double diff2 = z.b - w.b;
	long double diff3 = z.c - w.c;
	long double diff4 = z.d - w.d;

	// Scale by max component to avoid overflow/underflow
	long double max_diff = std::max(std::abs(diff1), std::max(std::abs(diff2), std::max(std::abs(diff3), std::abs(diff4))));
	
	if (max_diff == 0.0L) {
		return 0.0;
	}
	
	long double scaled1 = diff1 / max_diff;
	long double scaled2 = diff2 / max_diff;
	long double scaled3 = diff3 / max_diff;
	long double scaled4 = diff4 / max_diff;
	
	return max_diff * std::sqrt(scaled1*scaled1 + scaled2*scaled2 + scaled3*scaled3 + scaled4*scaled4);
}


int main(int argc, char const *argv[]) {

	auto ctx = prec::make_context("complex", argc, argv);
	ctx.output->settings.outputFiles = { "test/prec/prec_complex.csv" };
	random::random_source rnd = ctx.random->get_rnd();

	auto eq_opt = prec::equation_options<complex<>>(
		ctx.settings.defaultTolerance,
		distance_complex
	);

	auto phasor_opt = prec::equation_options<phasor<>>(
		ctx.settings.defaultTolerance,
		distance_phasor
	);

	auto quat_opt = prec::equation_options<quat<>>(
		ctx.settings.defaultTolerance,
		distance_quat
	);


	// complex.h

	// Complex number construction and basic properties
	{
		complex<> z = rand_complex(rnd);
		ctx.equals("complex::Re()", z.Re(), z.a);
		ctx.equals("complex::Im()", z.Im(), z.b);
	}

	// Complex number from real
	{
		complex<> z (rnd.uniform(-MAX, +MAX));
		ctx.equals("complex(real)", z, complex<>(z.a, 0), eq_opt);
	}

	// Zero complex number
	{
		complex<> z;
		ctx.equals("complex()", z, complex<>(0, 0), eq_opt);
	}

	// Complex conjugate
	{
		complex<> z = rand_complex(rnd);
		complex<> conj = z.conjugate();

		ctx.equals("complex::conjugate()", conj, complex<>(z.a, -z.b), eq_opt);
	}

	// Complex norm
	{
		complex<> z = rand_complex(rnd, 1E+03);
		ctx.equals("complex::sqr_norm() = complex::norm()**2", z.sqr_norm(), square(z.norm()));
	}

	// Complex argument (phase angle)
	{
		complex<> z1 (rnd.uniform(0, MAX), 0.0);
		ctx.equals("complex::arg()", z1.arg(), 0.0);

		complex<> z2 (0.0, rnd.uniform(0, MAX));
		ctx.equals("complex::arg()", z2.arg(), PI / 2.0);

		complex<> z3 (rnd.uniform(-MAX, 0), 0.0);
		ctx.equals("complex::arg()", th::abs(z3.arg()), PI);
	}

	// Complex inverse
	{
		complex<> z = rand_complex(rnd);
		complex<> inv = z.inverse();
		complex<> prod = z * inv;
		ctx.equals("complex::inverse()", prod.Re(), 1.0);
		ctx.equals("complex::inverse()", prod.Im(), 0.0);
	}

	// Complex addition
	{
		complex<> z1 = rand_complex(rnd);
		complex<> z2 = rand_complex(rnd);
		complex<> sum = z1 + z2;

		ctx.equals("complex::operator+()", sum.Re(), z1.a + z2.a);
		ctx.equals("complex::operator+()", sum.Im(), z1.b + z2.b);
	}

	// Complex subtraction
	{
		complex<> z1 = rand_complex(rnd);
		complex<> z2 = rand_complex(rnd);
		complex<> sub = z1 - z2;

		ctx.equals("complex::operator-()", sub.Re(), z1.a - z2.a);
		ctx.equals("complex::operator-()", sub.Im(), z1.b - z2.b);
	}

	// Complex multiplication
	{
		complex<> z1 (2.0, 3.0);
		complex<> z2 (4.0, 5.0);
		complex<> prod = z1 * z2;

		ctx.equals("complex::operator*()", prod.Re(), -7.0);
		ctx.equals("complex::operator*()", prod.Im(), 22.0);
	}

	// Complex division
	{
		complex<> z1 = rand_complex(rnd);
		complex<> z2 = rand_complex(rnd);
		complex<> quot = z1 / z2;

		ctx.equals("complex::operator/()", quot * z2, z1, eq_opt);
	}

	// Complex scalar multiplication
	{
		complex<> z = rand_complex(rnd);
		real scal = rnd.uniform(-MAX, +MAX);
		complex<> scaled = z * scal;

		ctx.equals("complex::operator*(real)", scaled, complex<>(z.a * scal, z.b * scal), eq_opt);
	}

	// Complex scalar division
	{
		complex<> z(6.0, 8.0);
		complex<> scaled = z / 2.0;
		ctx.equals("complex::operator/(real)", scaled, complex<>(3.0, 4.0), eq_opt);
	}

	// Complex negation
	{
		complex<> z(3.0, 4.0);
		complex<> neg = -z;
		ctx.equals("-complex", neg, complex<>(-3.0, -4.0), eq_opt);
	}

	// Complex equality
	{
		complex<> z1(3.0, 4.0);
		complex<> z2(3.0, 4.0);
		complex<> z3(5.0, 6.0);
		ctx.equals("complex::operator==()", z1 == z2, true);
		ctx.equals("complex::operator==()", z1 != z3, true);
	}

	// Complex assignment from array
	{
		complex<> z;
		std::array<real, 2> arr = {7.0, 8.0};
		z = arr;
		ctx.equals("complex({})", z, complex<>(7.0, 8.0), eq_opt);
	}


	// complex_analysis.h

	// Complex square
	{
		complex<> z = rand_complex(rnd);
		complex<> sqr = square(z);
		
		ctx.equals("square(complex)", sqr, z * z, eq_opt);
	}

	// Complex cube
	{
		complex<> z = rand_complex(rnd, 1E+03);
		complex<> cb = cube(z);

		ctx.equals("cube(complex)", cb, z * z * z, eq_opt);
	}

	// Complex exponential
	{
		complex<> z = complex<>(rnd.uniform(-10, 10), rnd.uniform(-10, 10));
		complex<> exp_z = th::exp(z);

		ctx.equals("exp(complex)", exp_z, complex<>(std::cos(z.b), std::sin(z.b)) * std::exp(z.a), eq_opt);
	}

	// Complex exponential at origin
	{
		complex<> z(0.0, 0.0);
		complex<> exp_z = exp(z);

		ctx.equals("exp(complex)", exp_z, complex<>(1, 0), eq_opt);
	}

	// Complex sine
	{
		real x = rnd.uniform(-MAX, +MAX);
		ctx.equals("sin(complex)", th::sin(complex<>(x, 0)), complex<>(std::sin(x), 0), eq_opt);
	}

	// Complex cosine
	{
		real x = rnd.uniform(-MAX, +MAX);
		ctx.equals("cos(complex)", th::cos(complex<>(x, 0)), complex<>(std::cos(x), 0), eq_opt);
	}

	// Complex square root identity
	{
		complex<> z = rand_complex(rnd);
		complex<> sq = sqrt(z);
		complex<> squared = sq * sq;

		ctx.equals("sqrt(complex)", squared, z, eq_opt);
	}

	// Complex logarithm
	{
		complex<> z (rnd.uniform(0, MAX), rnd.uniform(0, MAX));
		complex<> w (rnd.uniform(0, MAX), rnd.uniform(0, MAX));

		ctx.equals("ln(z*w) = ln(z)+ln(w)", th::ln(z * w), th::ln(z) + th::ln(w), eq_opt);
	}

	// Complex power
	{
		complex<> z = rand_complex(rnd);
		
		ctx.equals("pow(complex, 3)", th::pow(z, 3), z * z * z, eq_opt);
		ctx.equals("pow(complex, -1)", th::pow(z, -1), z.inverse(), eq_opt);
	}

	// Complex conjugate function
	{
		complex<> z = rand_complex(rnd, 1E+04);
		complex<> conj = conjugate(z);

		ctx.equals("z*z.conjugate() == z.sqr_norm()", z * conj, complex<>(z.sqr_norm()), eq_opt);
	}

	// Complex absolute value
	{
		complex<> z = rand_complex(rnd);
		real mag = abs(z);

		ctx.equals("abs(complex) == complex::norm()", mag, z.norm());
	}


	// phasor.h

	// Phasor construction from modulus and phase
	{
		phasor<> p = rand_phasor(rnd);
		ctx.equals("phasor(real, real)", p, phasor<>(p.modulus, p.phase), phasor_opt);
	}

	// Phasor from complex number
	{
		complex<> z(3.0, 4.0);
		phasor<> p(z);
		ctx.equals("phasor(complex)", p.modulus, 5.0);
		ctx.equals("phasor(complex)", p.phase, th::atan2(4.0, 3.0));
	}

	// Phasor from real number
	{
		real r = rnd.gaussian(0, VARIANCE);
		phasor<> p (r);

		ctx.equals("phasor(real)", p, phasor<>(r), phasor_opt);
	}

	// Phasor from negative real number
	{
		phasor<> p(-5.0);
		ctx.equals("phasor(real)", p. modulus, 5.0);
		ctx.equals("phasor(real)", p.phase, PI);
	}

	// Phasor real part extraction
	{
		phasor<> p = rand_phasor(rnd);
		real re = p.Re();
		real expected = p.modulus * th::cos(p.phase);
		ctx.equals("phasor::Re()", re, expected);
	}

	// Phasor imaginary part extraction
	{
		phasor<> p = rand_phasor(rnd);
		real re = p.Im();
		real expected = p.modulus * th::sin(p.phase);
		ctx.equals("phasor::Im()", re, expected);
	}

	// Phasor norm
	{
		phasor<> p = rand_phasor(rnd);

		ctx.equals("phasor::norm()", p.norm(), p.modulus);
	}

	// Phasor conjugate
	{
		phasor<> p = rand_phasor(rnd);
		phasor<> conj = p.conjugate();

		ctx.equals("phasor::conjugate()", conj, phasor<>(p.modulus, -p.phase), phasor_opt);
	}

	// Phasor inverse
	{
		phasor<> p = rand_phasor(rnd);
		phasor<> inv = p.inverse();

		ctx.equals("phasor::inverse()", inv, phasor<>(1 / p.modulus, -p.phase), phasor_opt);
	}

	// Phasor multiplication
	{
		phasor<> p1 = rand_phasor(rnd);
		phasor<> p2 = rand_phasor(rnd);
		phasor<> prod = p1 * p2;

		ctx.equals("phasor::operator*()", prod, phasor<>(p1.modulus * p2.modulus, p1.phase + p2.phase), phasor_opt);
	}

	// Phasor division
	{
		phasor<> p1 = rand_phasor(rnd);
		phasor<> p2 = rand_phasor(rnd);
		phasor<> quot = p1 / p2;

		ctx.equals("phasor::operator*()", quot, phasor<>(p1.modulus / p2.modulus, p1.phase - p2.phase), phasor_opt);
	}

	// Phasor power
	{
		phasor<> p = rand_phasor(rnd);
		phasor<> p3 = p * p * p;

		ctx.equals("phasor::operator*()", p3, phasor<>(cube(p.modulus), 3 * p.phase), phasor_opt);
	}

	// Phasor to complex and back
	{
		complex<> z = rand_complex(rnd, 1E+02);
		phasor<> p(z);

		complex<> z_back = complex<>(p.Re(), p.Im());
		ctx.equals("complex(phasor(complex))", z_back, z, eq_opt);
	}


	// quat.h

	// Quaternion construction
	{
		quat<> q (1.0, 2.0, 3.0, 4.0);

		ctx.equals("quat(real, ...)", q.a, 1.0);
		ctx.equals("quat(real, ...)", q.b, 2.0);
		ctx.equals("quat(real, ...)", q.c, 3.0);
		ctx.equals("quat(real, ...)", q.d, 4.0);
	}

	// Quaternion from real
	{
		real r = rnd.gaussian(0, VARIANCE);
		quat<> q (r);

		ctx.equals("quat(real)", q, quat<>(r, 0, 0, 0), quat_opt);
	}

	// Zero quaternion
	{
		quat<> q;

		ctx.equals("quat()", q, quat<>(0.0), quat_opt);
	}

	// Quaternion real and imaginary parts
	{
		quat<> q = rand_quat(rnd);

		ctx.equals("quat.Re()", q.Re(), q.a);
		ctx.equals("quat.Im1()", q.Im1(), q.b);
		ctx.equals("quat.Im2()", q.Im2(), q.c);
		ctx.equals("quat.Im3()", q.Im3(), q.d);
	}

	// Quaternion norm
	{
		quat<> q = rand_quat(rnd, 1E+04);
		ctx.equals("quat::sqr_norm() == quat::norm()^2", q.norm() * q.norm(), q.sqr_norm());
	}

	// Quaternion conjugate
	{
		quat<> q = rand_quat(rnd);
		quat<> conj = q.conjugate();

		ctx.equals("quat::conjugate()", conj, quat<>(q.a, -q.b, -q.c, -q.d), quat_opt);
	}

	// Quaternion normalized
	/*{
		quat<> q(1.0, 0.0, 0.0, 0.0);
		quat<> normalized = q.normalized();
		ctx.equals("quat(1,0,0,0).normalized().norm()", normalized.norm(), 1.0);
	}*/

	// Quaternion addition
	{
		quat<> q1 = rand_quat(rnd);
		quat<> q2 = rand_quat(rnd);
		quat<> sum = q1 + q2;

		ctx.equals("quat::operator+()", sum, quat<>(q1.a + q2.a, q1.b + q2.b, q1.c + q2.c, q1.d + q2.d), quat_opt);
	}

	// Quaternion subtraction
	{
		quat<> q1 = rand_quat(rnd);
		quat<> q2 = rand_quat(rnd);
		quat<> diff = q1 - q2;
		
		ctx.equals("quat::operator-()", diff, quat<>(q1.a - q2.a, q1.b - q2.b, q1.c - q2.c, q1.d - q2.d));
	}

	// Quaternion scalar multiplication
	{
		real r = rnd.gaussian(0, VARIANCE);
		quat<> q = rand_quat(rnd);
		quat<> scaled = q * r;

		ctx.equals("quat::operator*(real)", scaled, quat<>(r * q.a, r * q.b, r * q.c, r * q.d), quat_opt);
	}

	// Quaternion scalar division
	{
		real r = rnd.gaussian(0, VARIANCE);
		quat<> q = rand_quat(rnd);
		quat<> scaled = q / r;

		ctx.equals("quat::operator/(real)", scaled, quat<>(q.a / r, q.b / r, q.c / r, q.d / r), quat_opt);
	}

	// Quaternion multiplication (non-commutative)
	{
		quat<> q1(0.0, 1.0, 0.0, 0.0); // i
		quat<> q2(0.0, 0.0, 1.0, 0.0); // j
		quat<> prod = q1 * q2; // i*j = k

		ctx.equals("i*j = k", prod, quat<>(0, 0, 0, 1), quat_opt);
	}

	// Quaternion multiplication property (q * conjugate(q) = norm^2)
	{
		quat<> q = rand_quat(rnd, 1E+03);

		quat<> conj = q.conjugate();
		quat<> prod = q * conj;
		real expected = q.sqr_norm();

		ctx.equals("q*q.conjugate() == q.norm()^2", prod, quat<>(expected), quat_opt);
	}

	// Quaternion inverse
	{
		quat<> q = rand_quat(rnd);
		quat<> inv = q.inverse();
		quat<> product = q * inv;

		ctx.equals("q*q.inverse() == 1", product, quat<>(1.0), quat_opt);
	}

	// Quaternion negation
	{
		quat<> q = rand_quat(rnd);
		quat<> neg = -q;
		ctx.equals("-quat", neg, quat<>(-q.a, -q.b, -q.c, -q.d), quat_opt);
	}

	// Quaternion equality
	/*{
		quat<> q1(1.0, 2.0, 3.0, 4.0);
		quat<> q2(1.0, 2.0, 3.0, 4.0);
		quat<> q3(5.0, 6.0, 7.0, 8.0);
		ctx.equals("quat == quat", q1 == q2, true);
		ctx.equals("quat != quat", q1 != q3, true);
	}*/

	// Quaternion assignment from array
	{
		quat<> q;
		std::array<real, 4> v = {1.0, 2.0, 3.0, 4.0};
		q = v;
		ctx.equals("quat(std::array)", q, quat<>(v[0], v[1], v[2], v[3]), quat_opt);
	}

	// Quaternion division
	{
		quat<> q1 = rand_quat(rnd);
		quat<> q2 = rand_quat(rnd);
		quat<> quot = q1 / q2;

		ctx.equals("quat::operator/()", quot * q2, q1, quat_opt);
	}

	// Quaternion j*i = -k
	{
		quat<> j(0.0, 0.0, 1.0, 0.0);
		quat<> i(0.0, 1.0, 0.0, 0.0);
		quat<> prod = j * i; // j*i = -k

		ctx.equals("j*i = -k", prod, quat<>(0, 0, 0, -1), quat_opt);
	}


	// complex_types.h

	// Complex type trait
	{
		ctx.equals("is_complex_type<complex<real>>", is_complex_type<complex<real>>::value, true);
		ctx.equals("is_complex_type<complex<int>>", is_complex_type<complex<int>>::value, true);
		ctx.equals("is_complex_type<std::vector<real>>", is_complex_type<std::vector<real>>::value, false);
		ctx.equals("is_complex_type<quat<real>>", is_complex_type<quat<real>>::value, false);
		ctx.equals("is_complex_type<real>", is_complex_type<real>::value, false);
	}

	// Bi-complex (complex of complex)
	{
		bicomplex<> bz(complex<>(1.0, 2.0), complex<>(3.0, 4.0));
		ctx.equals("bicomplex::Re().Re()", bz.a.a, 1.0);
		ctx.equals("bicomplex::Im().Re()", bz.b.a, 3.0);
	}


	// Estimators

	// Complex exponential
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(0.0, 2.0 * PI),
			prec::estimator::quadrature1D()
		);

		ctx.homogeneous(
			"th::exp(complex)",
			[](real x) { return (exp(complex<>(0.0, x)) - complex<>(th::cos(x), th::sin(x))).norm(); },
			opt
		);
	}

	{
		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(-10, 10), prec::interval(-MAX, +MAX) },
			prec::estimator::montecarlo<real>(ctx.random, 2)
		);

		ctx.homogeneous(
			"th::exp(complex)",
			[](std::vector<real> v) {

				auto expected = std::exp(v[0]) * complex<>(std::cos(v[1]), std::sin(v[1]));
				auto eval = th::exp(complex<>(v[0], v[1]));

				return (eval - expected).norm();
			},
			opt
		);
	}

	{
		auto opt = prec::estimate_options<real, std::vector<real>>(
			{ prec::interval(-1E+07, 1E+07), prec::interval(-1E+07, +1E+07) },
			prec::estimator::montecarlo<real>(ctx.random, 2)
		);

		ctx.homogeneous(
			"th::sqrt(complex)",
			[](std::vector<real> v) {

				auto z = complex<>(v[0], v[1]);
				auto sqrt_z = th::sqrt(z);

				return (z - sqrt_z * sqrt_z).norm();
			},
			opt
		);
	}

	// Complex sine
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(-PI, PI),
			prec::estimator::quadrature1D()
		);

		ctx.homogeneous(
			"th::sin(complex)",
			[](real x) { return (th::sin(complex<>(x, 0.0)) - complex<>(std::sin(x), 0.0)).norm(); },
			opt
		);
	}

	// Complex cosine
	{
		auto opt = prec::estimate_options<real, real>(
			prec::interval(-PI, PI),
			prec::estimator::quadrature1D()
		);

		ctx.homogeneous(
			"th::cos(complex)",
			[](real x) { return (th::cos(complex<>(x, 0.0)) - complex<>(std::cos(x), 0.0)).norm(); },
			opt
		);
	}

}
