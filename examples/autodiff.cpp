
///
/// @file autodiff.cpp Automatic differentiation example.
/// This example illustrates how to use multivariate
/// automatic differentiation using the autodiff module
/// and the multidual<N> class.
/// You can compile it using 'make autodiff'
///

#include <iostream>
#include "theoretica.h"
using namespace th;


// A function from R^2 to R
template<typename NumType>
NumType f(vec<NumType, 2> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	// The th namespace is needed if <cmath>
	// functions are global, to avoid ambiguity
	return th::sqrt(x * y) * th::tan(x);
}


// A function from R^2 to R^2
autodiff::dvec2 g(autodiff::dvec2 v) {

	const auto x = v[0];
	const auto y = v[1];

	return {
		th::sqrt(x * y),
		x / y
	};
}


int main() {

	vec2 v = {1, 2};

	// Compute common differential operators on f(x, y)

	// You can call the function as usual
	io::println("f(v) =", f(v));


	// Use dreal2 (vector of dual2 numbers) for first-order autodiff
	auto df = f<autodiff::dreal2>;
	io::println("grad(f) =", autodiff::gradient(df, v));
	io::println("div(f) =", autodiff::divergence(df, v));


	// Use dual2 numbers for second-order autodiff
	auto d2f = f<dual2>;
	io::println("laplacian(f) =", autodiff::laplacian(d2f, v));


	// Compute the Jacobian matrix of g(x, y)
	io::println("jacobian(g):");
	io::println(autodiff::jacobian(g, v));

	return 0;
}
