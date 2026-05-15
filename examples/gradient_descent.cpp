
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

    vec2 x = multi_minimize(f<autodiff::dreal2>, vec2({1, 1}));

    // When using fixed-size containers, you may need to
    // specify the size by template argument.

    io::println("min at x =", x);
    io::println("f(x) =", f(x));
}
