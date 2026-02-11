
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;
using namespace autodiff;


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


dreal p(dvec v) {
    return v * v;
}


dreal H(dvec v) {
    return v * v + 1000;
}


dvec V(dvec v) {
	return {
		1 / (v * v), 1 / (v * v), 1 / (v * v)
	};
}


dreal d1(dvec v) {
    return v[0] - 2 * v[1] + v[2];
}


int main(int argc, char const *argv[]) {


	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);


	auto ctx = prec::make_context("autodiff");
	ctx.settings.outputFiles = { "test/prec/prec_autodiff.h" };
	ctx.settings.defaultIterations = 1000;

	// autodiff.h

	// dual.h

	// dual2.h

	// multidual.h
}
