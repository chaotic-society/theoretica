# Theoretica
<!-- Home -->
<!-- ======== -->

![GitHub last commit](https://img.shields.io/github/last-commit/chaotic-society/theoretica) ![GitHub Workflow Status](https://img.shields.io/github/actions/workflow/status/chaotic-society/theoretica/test-windows.yml) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/0f4ae5dc6e1140ad855a3d6325d44b35)](https://app.codacy.com/gh/chaotic-society/theoretica/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=chaotic-society/theoretica&amp;utm_campaign=Badge_Grade)  [![Documentation](https://img.shields.io/badge/Doxygen-docs-blue?style=flat&cacheSeconds=https%3A%2F%2Fchaotic-society.github.io%2Ftheoretica%2F&link=https%3A%2F%2Fchaotic-society.github.io%2Ftheoretica%2F)](https://chaotic-society.github.io/theoretica)  [![License](https://img.shields.io/github/license/chaotic-society/theoretica)](https://choosealicense.com/licenses/lgpl-3.0/)

> A C++ math library for scientific computing with a simple and elegant interface.

Theoretica provides methods for **scientific computing**, statistical analysis of experimental data and numerical approximations. The aim of the project is to give simple and immediate access to powerful algorithms for scientific and engineering software and graphical applications. The library is tested using [Chebyshev](https://github.com/chaotic-society/chebyshev), a unit testing framework specifically developed for scientific and numerical software.

## Table of Contents
- [Features](#features)
- [Setup](#setup)
- [Examples](#examples)
- [Documentation](#documentation)
- [Contributing](#contributing)
- [Workflow](#workflow)
- [License](#license)

---

## Features

Some features of the library include:

- Cross-platform with x86 and OpenMP enhancements
- Common real and complex functions
- Numerical Linear Algebra
- Complex numbers (algebraic & exponential form) and **quaternions**
- Numerical calculus, from integrals to ODEs
- Multivariate **Automatic Differentiation** with differential operators
- Descriptive and inferential statistics and **Monte Carlo** methods
- Bezier curves and spline interpolation

Theoretica is constantly developed and improved with new ideas!

## Setup

Theoretica is a header-only library and has **no dependencies**, so you can include it in your projects straight-away! You can build all tests and example programs by running `make all` in the main folder, ensuring that it works on your machine. When using the library, you can include single headers or use `theoretica.h` which includes all library headers (you can define `THEORETICA_INCLUDE_BASE` to make it include only fundamental headers).

### Quickstart

You can compile this simple code to check your setup:

```cpp
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
```

## Examples

The `examples` folder contains simple programs that showcase usage of the library:

- [Statistics](https://github.com/chaotic-society/theoretica/blob/master/examples/statistics.cpp)
- [Chaotic Attractors](https://github.com/chaotic-society/theoretica/blob/master/examples/attractor.cpp)
- [Automatic differentiation](https://github.com/chaotic-society/theoretica/blob/master/examples/autodiff.cpp)
- [Sampling distributions](https://github.com/chaotic-society/theoretica/blob/master/examples/sampling.cpp)
- [Monte Carlo integration](https://github.com/chaotic-society/theoretica/blob/master/examples/montecarlo_integral.cpp)
- [Fitting data](https://github.com/chaotic-society/theoretica/blob/master/examples/logfit.cpp)
- [Automatic error propagation](https://github.com/chaotic-society/theoretica/blob/master/examples/error_propagation.cpp)
- [Histogram usage](https://github.com/chaotic-society/theoretica/blob/master/examples/histogram.cpp)


## Documentation

The documentation for the project is available at this [link](https://chaotic-society.github.io/theoretica). The documentation is written using Doxygen syntax alongside the source code and the website is automatically updated on each commit. The HTML documentation is also available for download in the `gh-pages` branch. The bibliography used during research for the library is listed in [BIBLIOGRAPHY.md](https://github.com/chaotic-society/theoretica/blob/master/BIBLIOGRAPHY.md). You may learn more about the design choices behind the library reading the [Software Specification](https://github.com/chaotic-society/Theoretica-Lab/blob/main/specification/Theoretica_Software_Structure_Specification.pdf).

## Contributing

Contributions are welcome and appreciated! Have a look at the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/CONTRIBUTING.md) to learn more about how you can help. Contributions include writing code and documentation, testing and researching algorithms.

## Workflow

[![Test on Linux](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml) [![Test on Windows](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml) [![Test on MacOS](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml)

Theoretica uses automated workflows for recurring tasks. On each commit to `master`, tests are run on Linux, Windows and MacOS, benchmarks are launched and documentation is built and deployed to the online website. This ensures that the library works correctly and the documentation is always up-to-date.

## License

The project is currently under the [GNU Lesser General Public License 3.0](https://github.com/chaotic-society/theoretica/blob/master/LICENSE). You may learn more about it [here](https://choosealicense.com/licenses/lgpl-3.0/).
