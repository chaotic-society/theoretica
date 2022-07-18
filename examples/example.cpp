
///
/// @file example.cpp Example file
///

#include "theoretica.h"

using namespace th;

int main() {
 
    vec3 v = {1, 2, 3};
    mat3 A = mat3::identity();
    vec3 w = A * v;
 
    return 0;
}
