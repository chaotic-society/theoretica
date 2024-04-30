
///
/// @file example.cpp Example file.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
using namespace th;

#include "algebra/sparse_vec.h"


int main() {
    
    sparse_vec<> v1 = sparse_vec<>({
        vec2({1, PI}), vec2({9, PHI})
    });

    sparse_vec<> v2 = sparse_vec<>({
        vec2({1, 3.0})
    });

    std::cout << v1 << std::endl;
    std::cout << v2 << std::endl;
    std::cout << v1 * v2 << std::endl;
    
    char c;
    std::cin >> c;
    return 0;
}
