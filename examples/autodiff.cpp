
///
/// @file autodiff.cpp Automatic differentiation example.
/// This example may be compiled using 'make autodiff'
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
template<typename NumType>
vec<NumType, 2> g(vec<NumType, 2> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	return {th::sqrt(x * y), x / y};
}


int main() {

	vec2 v = {1, 2};

	// Compute common differential operators on f(x, y)

	// You can call the function as usual
	std::cout << "f(v) = " << f(v) << std::endl;

	// And also automatically compute differential operators
	std::cout << "grad(f) = " << gradient(f, v) << std::endl;
	std::cout << "div(f) = " << divergence(f, v) << "\n" << std::endl;
	std::cout << "laplacian(f) = " << laplacian(f, v) << "\n" << std::endl;

	// When you apply differential operators, the function is cast to
	// accept dual, multidual or dual2 arguments and is then evaluated
	// at the given point.

	// You can also apply the operators without specifying the point,
	// in which case the function returns a lambda function that
	// uses the same method to compute the operators when called.

	// Compute the Jacobian matrix of g(x, y)
	// Note that you may need to specify the input and output size
	// for template deduction, unlike the other functions.
	std::cout << "jacobian(g):\n" << jacobian<2, 2>(g, v) << std::endl;

	return 0;
}
