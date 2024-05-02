
///
/// @file example.cpp Example file.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
using namespace th;

#include "algebra/sparse_vec.h"


int main() {
    
    sparse_vec<> v1 = {
        vec2({1, PI}), vec2({9, PHI})
    };

    sparse_vec<> v2 = {
        vec2({1, -PI})
    };

    std::cout << v1 << std::endl;
    std::cout << v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    std::cout << v1 + v2 << std::endl;

    auto v3 = v1 + v2;
    v3.trim();

    std::cout << v3 << std::endl;

    return 0;
}
