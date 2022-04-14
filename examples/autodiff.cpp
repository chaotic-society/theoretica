
///
/// @file autodiff.cpp Automatic differentiation example
///

#include "uroboro.h"
using namespace umath;


// A function from R^2 to R
template<typename NumType>
NumType f(vec<2, NumType> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	return umath::sqrt(x * y) * umath::tan(x);
}


// A function from R^2 to R^2
template<typename NumType>
vec<2, NumType> g(vec<2, NumType> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	return {umath::sqrt(x * y), x / y};
}


int main() {

	vec2 v = {1, 2};

	// Compute common differential operators on f(x, y)
	std::cout << "f(v) = " << f(v) << std::endl;
	std::cout << "grad(f) = " << gradient(f, v) << std::endl;
	std::cout << "div(f) = " << divergence(f, v) << "\n" << std::endl;

	// Compute the Jacobian matrix of g(x, y)
	// Note that you may need to specify the input and output size
	// for template deduction, unlike the othe functions.
	std::cout << "jacobian(g):\n" << jacobian<2, 2>(g, v) << std::endl;

	return 0;
}
