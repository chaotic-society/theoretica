
///
/// @file multivariate_minimization.cpp Gradient descent minimization example.
/// This example may be compiled using 'make multivariate_minimization'
///

#include <iostream>
#include "theoretica.h"
using namespace th;


// The function to minimize.
// Should be declared either as a template to allow
// the algorithms to automatically compute the derivatives
// using autodiff.
template<typename NumType>
NumType f(vec<2, NumType> v) {

    const NumType x = v[0];
    const NumType y = v[1];

    return -th::exp(-square(x - 3) - square(y - 2));
}


int main() {

    // Use the best available algorithm to find
    // a minimum of the function starting from
    // a guess.
    vec2 x = minimize<2>(f, {1, 1});

    std::cout << "min at x = " << x << std::endl;
    std::cout << "f(x) = " << f(x) << std::endl;
 
    return 0;
}
