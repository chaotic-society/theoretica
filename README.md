<img alt="license badge" src="https://img.shields.io/github/license/mattiaisgro/uroboro"> <img alt="last commit badge" src="https://img.shields.io/github/last-commit/mattiaisgro/uroboro"> <img alt="code size in bytes badge" src="https://img.shields.io/github/languages/code-size/mattiaisgro/uroboro"> <img alt="issues badge" src="https://img.shields.io/github/issues/mattiaisgro/uroboro">

# Uroboro Math Library
A header-only, **C++ mathematical library**, focused on **graphical** and **physical** applications. Includes **real and complex analysis** functions with x86 Assembly enhancements, **linear algebra** operations, **quaternions** and **statistical** functions. Many other features are also supported, see [Functionalities](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Functionalities) for an overview of the functionalities provided by the library.

## Dependencies
The library has no dependencies. Only the C++ Standard Library is needed to use it.

## Functionalities
This is an overview of the library's functionalities. For a more detailed list see [FEATURES.md](https://github.com/mattiaisgro/uroboro/blob/master/FEATURES.md)
- Real analysis (`real_analysis.h`)
- Linear algebra (vector and matrix operations) (`algebra/vec.h`, `algebra/mat.h`)
- Complex (in algebraic and exponential form) and quaternion numbers (`complex/complex.h`, `complex/quat.h`, `complex/phasor.h`)
- Complex analysis (`complex/complex_analysis.h`)
- Dual numbers and Automatic Differentiation (`other/dual.h`, `other/dual_functions.h`)
- Statistical functions (including Least Squares Linearization) (`statistics/statistics.h`)
- Probability distribution functions (`statistics/distributions.h`)
- Approximation of roots and extrema of real functions (`approx.h`)
- Derivative approximation (`calculus/derivation.h`)
- Integral approximation, including Runge-Kutta of 4th order (`calculus/integration.h`)
- Polynomial operations, including derivation and integration (`polynomial.h`)
- Interpolation between vector data, including Bezier curves (`other/interpolation.h`)

## Usage
The library is header-only, so it is only needed to include the proper header files in your program to use it. To simplify the usage of the library, the `uroboro.h` file automatically includes common headers for real and complex analysis and linear algebra.
- Make sure to define `UROBORO_X86` if  you are using the x86 architecture, so that the library will use Assembly enhancements.
- To include **all header files** `UROBORO_INCLUDE_ALL` may be defined before including `uroboro.h`.
- All functions and classes are defined inside the `uroboro` namespace (and eventually other sub-namespaces). Another namespace alias is automatically defined as `umath` (same as using `uroboro`). If you want to disable the `umath` alias, define `UROBORO_NO_NAMESPACE` before including `uroboro.h`

**Matrices** (`algebra/mat.h`) are implemented inside the `mat<N, K>` class (where N and K are the number of columns and rows respectively), while **vectors** (`algebra/vec.h`) in the `vec<N>` class (where `N` is the number of rows). **Quaternions** (`complex/quat.h`) and **complex** (`complex/complex.h`) numbers are implemented in the `quat` and `complex` classes. For example usage see [Examples](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Examples).

## Examples
Introductory examples can be found in [EXAMPLES.md](https://github.com/mattiaisgro/uroboro/blob/master/EXAMPLES.md), see  [example.cpp](https://github.com/mattiaisgro/uroboro/blob/master/src/example.cpp) for more example code.

### Quickstart
To check if you correctly set up the library, try to compile the following code in a source file:
```c
#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace uroboro;

int main() {
 
    vec3 v = {1, 2, 3};
    mat3 A = mat3::rotation_x_3x3(PI2);
 
    vec3 w = A * v;
 
    complex z = complex::rotor(PI4);
 
    quat q = {1, 2, 3, 4};
    q = q.inverse();
 
    return 0;
}
```

## Error handling
The library uses `errno` and `umath::MathException` (if it is enabled) to report errors. The behaviour of the library may be modified using the `UROBORO_THROW_EXCEPTIONS` and `UROBORO_ONLY_EXCEPTIONS`. See [Macros](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Macros) to learn more.

## Progress
Most functionalities have already been implemented, but test cases and documentation are in need of writing. Benchmarking different algorithms may be also used to enchance the library. Core functions which are currently under development are:
- Matrix determinant
- powf, exp, sin, cos, tan (only the **architecture independent** versions, they already work properly on x86 machines)
- asin, acos, atan, atan2 (Inverse trig functions)

See [Future development](https://github.com/mattiaisgro/uroboro/blob/master/README.md#future-development) for information about new features which might get implemented.

## Macros
These are the macros that can be defined to change the library's behaviour:
- **UROBORO_INCLUDE_ALL** - Including `uroboro.h` will include _all_ header files
- **UROBORO_BRANCHLESS** - Branchless implementations will be preferred when possible (and computationally convenient)
- **UROBORO_FORCE_BRANCHLESS** - Branchless implementations will always be used whenever possible
- **UROBORO_THROW_EXCEPTIONS** - Exceptions will be thrown and errno set on error (by default errno is set and NaN is returned)
- **UROBORO_ONLY_EXCEPTIONS** - Exceptions will be thrown on error (without modifying errno)
- **UROBORO_NO_NAMESPACE** - The `umath` namespace alias will not be defined
- **UROBORO_X86** - Assembly x86 implementations will be used whenever possible
- **UROBORO_FLOAT_PREC** - Floating point precision (`float`) will be used for the `real` type (by default `double` is used)
- **UROBORO_LONG_DOUBLE_PREC** - Long double precision (`long double`) will be used
- **UROBORO_ARBITRARY_PREC** - Arbitrary precision will be used (NOT implemented yet)
- See `constants.h` for more specific defines

## Future development
### In progress
These are features which are currently under development:
- Matrix determinant
- Column/row precedence independence for matrices (let the user decide whether to use column/row-first storage and colexicographical/lexicographical representation)

### Planned
These are features which will be developed soon:
- Romberg integral approximation
- Pseudorandom Number Generation
- LU decomposition of matrices
- Matrix eigenvalue finding
- Montecarlo integral approximation
- Benchmarking (either implemented or integrated)
- Make the library independent from the C++ Standard Library (at least partially)
- Metropolis pseudorandom number generation
- Noise generation

### Potential
These are interesting features which might or might not be implemented in the future.
- Sorting algorithms
- SSE implementations of vector and matrix operations

## Contributing
Pull requests are well accepted, open an issue or discussion for feature requests, suggestions or to point out problems within the library.

## License
This project is currently under the GNU Lesser General Public License 3.0




