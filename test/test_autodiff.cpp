
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


vec<> rand_vec(unsigned int n, interval I, PRNG& g) {

	vec<> v(n);
	for (unsigned int i = 0; i < n; ++i)
		v[i] = rand_uniform(I.a, I.b, g);

	return v;
}


// Check that a given operator applied to a certain function is equal to zero
template<typename Function>
prec::estimate_result test_operator(
	Function f, PRNG g, interval I, Real tol, unsigned int n, unsigned int size) {

	Real max = 0, sum = 0, sum2 = 0;

	for (size_t i = 0; i < n; ++i) {
		
		real res = f(rand_vec(size, I, g));

		max = th::max(max, th::abs(res));
		sum += th::abs(res);
		sum2 += square(res);
	}

	prec::estimate_result res;
	res.max_err = max;
	res.abs_err = sum;
	res.rms_err = th::sqrt(sum2) / n;
	res.mean_err = sum / n;
	res.rel_err = 0;

	if(res.max_err > tol)
		res.failed = true;

	return res;
}


dual2 h1(vec<dual2> v) {
	return ln(v[0] * v[0] + v[1] * v[1]);
}


dual2 h2(vec<dual2> v) {
	return exp(v[0]) * sin(v[1]);
}


d_real<> f(d_vec<> v) {
    return v * v;
}


d_real<> H(d_vec<> v) {
    return v * v + 1000;
}


d_vec<> V(d_vec<> v) {
	return {
		1 / (v * v), 1 / (v * v), 1 / (v * v)
	};
}


d_real<> d1(d_vec<> v) {

    return v[0] - 2 * v[1] + v[2];
}


int main(int argc, char const *argv[]) {


	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);


	prec::state.outputFolder = "test/";
	prec::state.defaultIterations = 1000;
	
	prec::setup("autodiff");

		// Compare the automatic derivative to the analytical derivative
		prec::estimate("dual::Dual()",
			[](real x) {
				dual d = dual(x, 1);
				return (th::cos(square(d)) / th::exp(-square(d)) / ln(1 / square(d))).Dual();
			},
			[](real x) {
				return (2 * th::exp(square(x)) * ((square(x) * ln(1 / square(x)) + 1)
							* th::cos(square(x)) - square(x) * ln(1 / square(x)) * th::sin(square(x))))
								/ (x * square(ln(1 / square(x))));
			}, {interval(0.001, 0.5), interval(-0.5, -0.001)}
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

		
		prec::estimate("dual::gradient (f)",
			[=](interval k, Real tol, unsigned int n) {
				return test_operator(
					[=](vec<> v) {
						return gradient(f, v) * (J * gradient(H, v));
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
