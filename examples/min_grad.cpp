
///
/// @file min_grad.cpp Gradient descent minimization example
///

#include "theoretica.h"
using namespace th;


template<typename NumType>
NumType f(vec<2, NumType> v) {

    const NumType x = v[0];
    const NumType y = v[1];

    return -exp(-square(x - 3) - square(y - 2));
}


int main() {

    vec2 x = minimize<2>(f, {1, 1});

    std::cout << "min at x = " << x << std::endl;
    std::cout << "f(x) = " << f(x) << std::endl;
 
    return 0;
}
