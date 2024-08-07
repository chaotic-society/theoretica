
///
/// @file gradient_descent.cpp Gradient descent minimization example.
/// This example may be compiled using 'make gradient_descent'
///

#include <iostream>
#include "theoretica.h"
using namespace th;


// The function to minimize.
// Should be declared as a template to allow
// the algorithms to automatically compute the derivatives
// using autodiff.
template<typename NumType>
NumType f(vec<NumType, 2> v) {

    const NumType x = v[0];
    const NumType y = v[1];

    return -th::exp(-square(x - 3) - square(y - 2));
}


int main() {

    // Use the best available algorithm to find
    // a minimum of the function starting from
    // a guess of x = 1 and y = 1.

    vec2 x = multi_minimize<2>(f, {1, 1});

    // When using fixed-size containers, you may need to
    // specify the size by template argument.

    std::cout << "min at x = " << x << std::endl;
    std::cout << "f(x) = " << f(x) << std::endl;
}
