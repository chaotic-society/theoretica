# Uroboro Mathematical Library

<img alt="" src="https://img.shields.io/github/license/mattiaisgro/uroboro"> <img alt="" src="https://img.shields.io/github/last-commit/mattiaisgro/uroboro"> <img alt="" src="https://img.shields.io/github/issues/mattiaisgro/uroboro">

A C++ numerical and automatic **mathematical library**, focused on **graphical** and **physical** applications. Includes **real** and **complex analysis** functions with x86 Assembly enhancements, **linear algebra** operations, **quaternions**, **statistical** functions and **numerical** and **automatic** methods for real functions. Many other features are also supported, see [Functionalities](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Functionalities) for an overview of the functionalities provided by the library.

## Dependencies
The library has no dependencies. Only the C++ Standard Library is needed to use it.

## Functionalities
This is an overview of the library's functionalities. For a more detailed list see [FEATURES.md](https://github.com/mattiaisgro/uroboro/blob/master/txt/FEATURES.md)
- **Real analysis**
- **Linear algebra** with vector and matrix operations
- **Complex** numbers (in algebraic and exponential form) and complex analysis functions
- **Quaternions**
- Dual numbers, **Multivariable Automatic Differentiation** and Differential Operators
- **Statistical** functions, including Least Squares Linearization
- Probability distribution functions
- Pseudorandom and Quasirandom number generation (LCG, Xoshiro256++, Splitmix64, Wyrand, Weyl)
- Random distribution sampling
- Approximation of **roots** and **extrema** of real functions
- Derivative approximation
- Integral approximation, including **Runge-Kutta** of 4th order and **Romberg integration**
- Polynomial operations and **polynomial interpolation** with Chebyshev nodes
- Spline curves, including generic Bezier curves

## Usage
The library is header-only, so it is only needed to include the proper header files in your program to use it. To simplify the usage of the library, the `uroboro.h` file automatically includes common headers for real and complex analysis and linear algebra.
- Make sure to define `UROBORO_X86` if  you are using the x86 architecture, so that the library will use Assembly enhancements.
- To include **all header files** `UROBORO_INCLUDE_ALL` may be defined before including `uroboro.h`.
- All functions and classes are defined inside the `uroboro` namespace (and eventually other sub-namespaces). Another namespace alias is automatically defined as `umath` (same as using `uroboro`). If you want to disable the `umath` alias, define `UROBORO_NO_NAMESPACE` before including `uroboro.h`

**Matrices** (`algebra/mat.h`) are implemented inside the `mat<N, K>` class (where N and K are the number of rows and columns respectively), while **vectors** (`algebra/vec.h`) in the `vec<N>` class (where `N` is the number of rows). **Quaternions** (`complex/quat.h`) and **complex** (`complex/complex.h`) numbers are implemented in the `quat` and `complex` classes. For example usage see [Examples](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Examples).

## Examples
Introductory examples can be found in [EXAMPLES.md](https://github.com/mattiaisgro/uroboro/blob/master/txt/EXAMPLES.md), see  [example.cpp](https://github.com/mattiaisgro/uroboro/blob/master/examples/example.cpp) for more example code.

### Quickstart
To check if you correctly set up the library, try to compile the following code in a source file:
```c
#define UROBORO_INCLUDE_ALL
#include "uroboro.h"

using namespace uroboro;

int main() {
 
    vec3 v = {1, 2, 3};
    mat3 A = mat3::rotation_3x3(PI4, {1, 0, 1});
    vec3 w = A * v;
 
    complex z = complex::rotor(PI4);
 
    quat q = {1, 2, 3, 4};
    q = q.inverse();
 
    return 0;
}
```

## Error handling
The library uses `errno` and `umath::MathException` (if it is enabled) to report errors. The behaviour of the library may be modified using the `UROBORO_THROW_EXCEPTIONS` and `UROBORO_ONLY_EXCEPTIONS`. See [Macros](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Macros) to learn more.

## Macros
These are the macros that can be defined to change the library's behaviour:
- **UROBORO_INCLUDE_ALL** - Including `uroboro.h` will include _all_ header files
- **UROBORO_BRANCHLESS** - Branchless implementations will be preferred when possible (and computationally convenient)
- **UROBORO_FORCE_BRANCHLESS** - Branchless implementations will always be used whenever possible
- **UROBORO_THROW_EXCEPTIONS** - Exceptions will be thrown and errno set on error (by default errno is set and NaN is returned)
- **UROBORO_ONLY_EXCEPTIONS** - Exceptions will be thrown on error (without modifying errno)
- **UROBORO_NO_NAMESPACE** - The `umath` namespace alias will **not** be defined
- **UROBORO_X86** - **Assembly x86** implementations will be used whenever possible
- **UROBORO_FLOAT_PREC** - Floating point precision (`float`) will be used for the `real` type (by default `double` is used)
- **UROBORO_LONG_DOUBLE_PREC** - Long double precision (`long double`) will be used
- **UROBORO_ARBITRARY_PREC** - Arbitrary precision will be used (NOT implemented yet)
- **UROBORO_NO_PRINT** - Do **not** compile `to_string()` and `operator<<()` methods for all classes (to avoid including `<string>`, `<sstream>` and `<ostream>`)
- **UROBORO_ROW_FIRST** - The `mat<N, K>` class will use row-first storage of matrix data instead of column-first.
- **UROBORO_MATRIX_LEXIC** - Lexicographical notation for matrices (column first access) will be used for matrix functions `at`, `get` and `set`. By default, matrix indices refer to row and column, in this order.
- See `constants.h` for more specific defines

## Contributing
Contributions are welcome and appreciated, please read the [Contributing Guide](https://github.com/chaotic-society/uroboro/blob/master/CONTRIBUTING.md) to learn about how to contribute to the project.

## License
This project is currently under the [GNU Lesser General Public License 3.0](https://github.com/chaotic-society/uroboro/blob/master/LICENSE).




