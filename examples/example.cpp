
///
/// @file example.cpp Example file.
/// This example may be compiled using 'make example'
///

#include <iostream>
#include <ctime>
#include <iomanip>

#include "theoretica.h"
#include "utility.h"
using namespace th;

autodiff::dreal f_multidual(autodiff::dvec v) {

    return {
        v[0] + v[1]
    };
}

dual f_dual(dual x) {
    return x;
}

dual2 f_dual2(dual2 x) {
    return x;
}

real f_real(real x) {
    return x;
}


complex<> f_complex(complex<> x) {
    return x;
}

struct st {};

real f_not(st s) {
    return 0.0;
}


int main() {


    TH_DEBUG(has_multidual_return<decltype(f_real)>());
    TH_DEBUG(has_multidual_return<decltype(f_not)>());
    TH_DEBUG(has_multidual_return<decltype(f_multidual)>());

    TH_DEBUG(autodiff::gradient(f_multidual, vec<real>({0, 0})));

}
