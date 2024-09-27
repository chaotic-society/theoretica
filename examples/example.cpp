
///
/// @file example.cpp Example program.
/// This example may be compiled using 'make example'
///

#include "theoretica.h"
using namespace th;

int main() {
    
    // Declare a 3D vector
    vec3 v = {1, 2, 3};

    // Create a 3x3 identity matrix
    mat3 A = mat3::identity();

    // Transform v by A
    vec3 w = A * v;
}
