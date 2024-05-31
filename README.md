# Theoretica

<!-- Home -->
<!-- ======== -->

![GitHub last commit](https://img.shields.io/github/last-commit/chaotic-society/theoretica) ![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/chaotic-society/theoretica/test-windows.yml) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/0f4ae5dc6e1140ad855a3d6325d44b35)](https://app.codacy.com/gh/chaotic-society/theoretica/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=chaotic-society/theoretica&amp;utm_campaign=Badge_Grade)  ![License](https://img.shields.io/github/license/chaotic-society/theoretica)

> A numerical and automatic math library for scientific research and graphical applications

Theoretica is a header-only mathematical library which provides algorithms for **systems simulation**, **statistical analysis** of experimental data and **numerical approximation**, using a **functional** oriented paradigm to mimic **mathematical notation** and formulas. The aim of the library is to provide _simple_ access to powerful algorithms while keeping an _elegant_ and _transparent_ interface, enabling the user to focus on the problem at hand.

### A short example
Given a Hamiltonian function H(q, p) and a function f(q,  p) defined on its phase space, you can compute its _exact_ time derivative at a position eta = (q, p) like this:
$$\frac{df}{dt} = \nabla f(\vec \eta) \cdot J \cdot \nabla H(\vec \eta)$$

Where J is the symplectic matrix. The equation can be translated into code as:
```java
mat<N, N> J = mat<N, N>::symplectic();
real df_dt = gradient(f, eta) * J * gradient(H, eta);
```

The library includes real and complex analysis functions optimized for the x86 architecture, linear algebra, quaternions, roots and extrema search, numerical approximation of derivatives, integrals and differential equations, as well as more advanced features like dual numbers for automatic differentiation, statistical functions including distribution sampling, pseudorandom and quasirandom number generation for Monte Carlo methods and simulations.

## Table of Contents
- [Key Features](#key-features)
- [Dependencies](#dependencies)
- [Setup](#setup)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Testing](#testing)
- [Branches](#branches)
- [Other information](#other-information)

---

## Key Features
This is an overview of the library's functionalities. For a more detailed list see [FEATURES.md](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/FEATURES.md)
- **Real and complex analysis**
- **Linear algebra** with common vector and matrix operations
- Complex numbers in algebraic and exponential form
- **Quaternions**
- Dual numbers, **Multivariable Automatic Differentiation** and Differential Operators
- Pseudorandom and Quasirandom number generation (LCG, Xoshiro256++, Splitmix64, Wyrand, Weyl)
- Statistical functions, including **Least Squares Linearization**
- Random distribution sampling and Monte Carlo
- Approximation of roots, extrema, derivatives and integrals of real functions
- Numerical integration of **Ordinary Differential Equations** (Euler, Heun, RK4, K38, multistep)
- Polynomial interpolation with Chebyshev nodes, Bezier curves and spline interpolation

## Dependencies
The library has no dependencies. Only the C++ Standard Library with C++11 capabilities is needed to use it. **You can include it in your project straight away!**

## Setup
You don't need anything other than your compiler to use the library. You can run `make all` in the root directory of the library to make sure it works.
Define **THEORETICA_INCLUDE_BASE** if you intend to use only basic functionalities (linear algebra, real functions, complex numbers), as by default `theoretica.h` includes all headers.
All library functions and objects are implemented in the `theoretica` namespace (`th` is a shorter namespace alias).

## Documentation
The documentation for the library is available [here](https://chaotic-society.github.io/theoretica/).
Introductory examples can be found in [EXAMPLES.md](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/EXAMPLES.md) and more advanced examples can be found inside the `examples/` folder.
The bibliography used for researching the algorithms used in the library is available [here](https://github.com/chaotic-society/theoretica/blob/master/docs/txt/BIBLIOGRAPHY.md).

### Quickstart
You can try to compile this simple code to check if you set up the library correctly:
```cpp
#include "theoretica.h"
using namespace th;

int main() {
 
    vec3 v = {1, 2, 3};
    mat3 A = mat3::identity();
    vec3 w = A * v;
 
    return 0;
}
```

### Other examples
- [Statistics](https://github.com/chaotic-society/theoretica/blob/master/examples/statistics.cpp)
- [Chaotic Attractors](https://github.com/chaotic-society/theoretica/blob/master/examples/attractor.cpp)
- [Automatic differentiation](https://github.com/chaotic-society/theoretica/blob/master/examples/autodiff.cpp)
- [Sampling distributions](https://github.com/chaotic-society/theoretica/blob/master/examples/sampling.cpp)
- [Monte Carlo integration](https://github.com/chaotic-society/theoretica/blob/master/examples/montecarlo_integral.cpp)
- [Fitting data](https://github.com/chaotic-society/theoretica/blob/master/examples/logfit.cpp)
- [Automatic error propagation](https://github.com/chaotic-society/theoretica/blob/master/examples/error_propagation.cpp)
- [Histogram usage](https://github.com/chaotic-society/theoretica/blob/master/examples/histogram.cpp)

## Contributing
Contributions are welcome and very appreciated! Make sure to read the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/CONTRIBUTING.md) to know more about how you can help. If you participate, you are expected to follow the [Code of Conduct](https://github.com/chaotic-society/theoretica/blob/master/CODE_OF_CONDUCT.md).

## Testing
[![Test on Linux](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml) [![Test on Windows](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml) [![Test on MacOS](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml)

The library uses the custom built [Chebyshev testing framework](https://github.com/chaotic-society/chebyshev) to **estimate the precision** of functions and test their **performance**. Tests are automatically run on Windows, Linux and MacOS on every commit to ensure stability. Test units are placed inside the `test` folder while benchmarks are placed inside the `benchmark` folder.

## Branches
- ### `main`
Where the current stable version of the library is developed.
- ### `physics`
The physics module implementing quantum mechanics features.
- ### `sparse`
Sparse matrix and vector algebra classes.
- ### `machine-learning`
Machine learning development branch, with cost functions, models and training algorithms.

## Other information

### License
This project is currently under the [GNU Lesser General Public License 3.0](https://github.com/chaotic-society/theoretica/blob/master/LICENSE).

### Macros
These are common macros that can be defined to change the library's behaviour:
| Macro | Description |
| ----- | ----------- |
|**THEORETICA_INCLUDE_BASE**|Including `theoretica.h` will only include base headers|
|**THEORETICA_THROW_EXCEPTIONS**|Exceptions will be thrown and errno set on error (by default errno is set and NaN is returned)|
|**THEORETICA_ONLY_EXCEPTIONS**|Exceptions will be thrown on error (without modifying errno)|
|**THEORETICA_X86**|**Assembly x86** implementations will be used whenever possible (automatically defined on most compilers by the library)|
|**THEORETICA_FLOAT_PREC**|Floating point precision (`float`) will be used for the `real` type (by default `double` is used)|
|**THEORETICA_LONG_DOUBLE_PREC**|Long double precision (`long double`) will be used|
|**THEORETICA_ROW_FIRST**|The `mat<N, K>` class will use row-first storage of matrix data instead of column-first.|

See `src/core/constants.h` for more specific macros.

### Error handling
The library uses `errno` and `th::math_exception` (if enabled) to report errors. The behaviour of the library may be modified using the `THEORETICA_THROW_EXCEPTIONS` and `THEORETICA_ONLY_EXCEPTIONS`. See [Macros](#Macros) to learn more.
