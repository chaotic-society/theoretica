# Theoretica

<img alt="" src="https://img.shields.io/github/last-commit/mattiaisgro/theoretica"> [![Test on Linux](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml) [![Test on Windows](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml) [![Test on MacOS](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/0f4ae5dc6e1140ad855a3d6325d44b35)](https://www.codacy.com/gh/chaotic-society/theoretica/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=chaotic-society/theoretica&amp;utm_campaign=Badge_Grade)  <img alt="" src="https://img.shields.io/github/license/mattiaisgro/theoretica">

A numerical and automatic mathematical library in C++ for scientific and graphical applications, using a mixed functional/object-oriented paradigm to **mimic elegant mathematical notation**.

### A short example
Given a Hamiltonian function **H(q, p)** and a function **f(q, p)** defined on its phase space, you can compute its _exact_ time derivative at a position _eta_ like this:
```cpp
real df_dt = gradient(f, eta) * mat<M, M>::symplectic() * gradient(H, eta);
```

The library includes basic functionalities like real and complex analysis functions enhanced for x86 architectures, vector and matrix operations (both row-major and column-major), quaternions and function roots and extrema search, as well as more advanced features like dual numbers for automatic differentiation, statistical functions including distribution sampling, pseudorandom and quasirandom number generation for Monte Carlo methods and simulations.

## Dependencies
The library has no dependencies. Only the C++ Standard Library is needed to use it. **You can include it in your project straight away!**

## Key Features
This is an overview of the library's functionalities. For a more detailed list see [FEATURES.md](https://github.com/mattiaisgro/theoretica/blob/master/txt/FEATURES.md)
- **Real and complex analysis**
- **Linear algebra** with common vector and matrix operations
- Complex numbers in algebraic and exponential form
- **Quaternions**
- Dual numbers, **Multivariable Automatic Differentiation** and Differential Operators
- Pseudorandom and Quasirandom number generation (LCG, Xoshiro256++, Splitmix64, Wyrand, Weyl)
- Statistical functions, including **Least Squares Linearization**
- Random distribution sampling
- Approximation of roots, extrema, derivatives and integrals of real functions
- Polynomial interpolation with Chebyshev nodes, Bezier curves and spline interpolation

## Setup
You don't need anything other than your compiler to use the library. You can run `make all` in the root directory of the library to make sure it works. 
- Define **THEORETICA_X86** if you are on an x86 architecture to enable Assembly instructions
- Define **THEORETICA_INCLUDE_ALL** if you intend to use most functionalities, as by default `theoretica.h` only includes base headers

All library functions are implemented in the `theoretica` namespace (`th` is a shorter namespace alias).

## Examples
Introductory examples can be found in [EXAMPLES.md](https://github.com/mattiaisgro/theoretica/blob/master/txt/EXAMPLES.md) and more advanced examples can be found inside the `examples/` folder.

### Quickstart
You can try to compile this simple code to see if the library works properly:
```cpp
#define THEORETICA_INCLUDE_ALL
#include "theoretica.h"

using namespace theoretica;

int main() {
 
    vec3 v = {1, 2, 3};
    mat3 A = mat3::identity();
    vec3 w = A * v;
 
    return 0;
}
```

## Contributing
Contributions are welcome and very appreciated! Make sure to read the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/CONTRIBUTING.md) to know more about how you can help. If you participate, you are expected to follow the [Code of Conduct](https://github.com/chaotic-society/theoretica/blob/master/CODE_OF_CONDUCT.md).

## License
This project is currently under the [GNU Lesser General Public License 3.0](https://github.com/chaotic-society/theoretica/blob/master/LICENSE).

## Macros
These are the macros that can be defined to change the library's behaviour:
- **THEORETICA_INCLUDE_ALL** - Including `theoretica.h` will include _all_ header files instead of base headers
- **THEORETICA_THROW_EXCEPTIONS** - Exceptions will be thrown and errno set on error (by default errno is set and NaN is returned)
- **THEORETICA_ONLY_EXCEPTIONS** - Exceptions will be thrown on error (without modifying errno)
- **THEORETICA_X86** - **Assembly x86** implementations will be used whenever possible
- **THEORETICA_FLOAT_PREC** - Floating point precision (`float`) will be used for the `real` type (by default `double` is used)
- **THEORETICA_LONG_DOUBLE_PREC** - Long double precision (`long double`) will be used
- **THEORETICA_ROW_FIRST** - The `mat<N, K>` class will use row-first storage of matrix data instead of column-first.
- **THEORETICA_MATRIX_LEXIC** - Lexicographical notation for matrices (column first access) will be used for matrix functions `at`, `get` and `set`. By default, matrix indices refer to row and column, in this order.
- See `constants.h` for more specific defines

## Error handling
The library uses `errno` and `th::MathException` (if it is enabled) to report errors. The behaviour of the library may be modified using the `THEORETICA_THROW_EXCEPTIONS` and `THEORETICA_ONLY_EXCEPTIONS`. See [Macros](#Macros) to learn more.
