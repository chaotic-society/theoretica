
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


real f(real x) {
	dual d = dual(x, 1);
	return (th::cos(square(d)) / th::exp(-square(d)) / ln(1 / square(d))).Dual();
}


real df(real x) {
	return (2 * th::exp(square(x)) * ((square(x) * ln(1 / square(x)) + 1)
		* th::cos(square(x)) - square(x) * ln(1 / square(x)) * th::sin(square(x))))
			/ (x * square(ln(1 / square(x))));
}


// Harmonic functions
dual2 h1(vec<dual2> v) {
	return ln(v[0] * v[0] + v[1] * v[1]);
}


dual2 h2(vec<dual2> v) {
	return exp(v[0]) * sin(v[1]);
}


dreal_t<> p(dvec_t<> v) {
    return v * v;
}


dreal_t<> H(dvec_t<> v) {
    return v * v + 1000;
}


dvec_t<> V(dvec_t<> v) {
	return {
		1 / (v * v), 1 / (v * v), 1 / (v * v)
	};
}


dreal_t<> d1(dvec_t<> v) {
    return v[0] - 2 * v[1] + v[2];
}


int main(int argc, char const *argv[]) {


	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);


	prec::state.outputFolder = "test/";
	prec::state.defaultIterations = 1000;
	
	prec::setup("autodiff");

		// Compare the automatic derivative to the analytical derivative
		prec::estimate(
			"dual::Dual()", f, df, {interval(0.001, 0.5), interval(-0.5, -0.001)}
		);


		prec::estimate("dual::laplacian (h1)",
			[](interval k, Real tol, unsigned int n) {
				return test_operator(
					laplacian(h1),
					PRNG::xoshiro(time(nullptr)),
					k, 1E-8, n, 2
				);
			}, interval(-100, 100)
		);


		prec::estimate("dual::laplacian (h2)",
			[](interval k, Real tol, unsigned int n) {
				return test_operator(
					laplacian(h2),
					PRNG::xoshiro(time(nullptr)),
					k, 1E-8, n, 2
				);
			}, interval(-100, 100)
		);


		// Test the gradient computation by evaluating
		// the time derivative of a constant of motion
		// for a Hamiltonian system
		mat<> J = mat4::symplectic();

		
		prec::estimate("dual::gradient (p)",
			[=](interval k, Real tol, unsigned int n) {
				return test_operator(
					[=](vec<> v) {
						return gradient(p, v) * (J * gradient(H, v));
					},
					PRNG::xoshiro(time(nullptr)),
					k, 1E-8, n, 4
				);
			}, interval(-100, 100)
		);


		// Test against an irrotational vector field
		prec::estimate("dual::curl (V)",
			[=](interval k, Real tol, unsigned int n) {
				return test_operator(
					[=](vec<> v) {
						return curl(V, v).sqr_norm();
					},
					PRNG::xoshiro(time(nullptr)),
					k, 1E-5, n, 3
				);
			}, interval(1, 100)
		);


		// Test against a divergence-free scalar field
		prec::estimate("dual::divergence (d1)",
			[](interval k, Real tol, unsigned int n) {
				return test_operator(
					divergence(d1),
					PRNG::xoshiro(time(nullptr)),
					k, 1E-8, n, 3
				);
			}, interval(-100, 100)
		);


	prec::terminate();
}
