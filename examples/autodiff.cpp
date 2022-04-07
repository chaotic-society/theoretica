
///
/// @file autodiff.cpp Automatic differentiation example
///

#include "uroboro.h"
using namespace umath;


template<typename NumType>
NumType f(vec<2, NumType> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	return umath::sqrt(x * y) * umath::tan(x);
}


int main() {

	vec2 v = {1, 1};

	std::cout << "f(v) = " << f(v) << std::endl;
	std::cout << "grad(f) = " << gradient(f, v) << std::endl;
	std::cout << "div(f) = " << divergence(f, v) << std::endl;

	return 0;
}
