# ⚛️ Theoretica
<!-- Home -->
<!-- ======== -->

![GitHub last commit](https://img.shields.io/github/last-commit/chaotic-society/theoretica) [![Codacy Badge](https://app.codacy.com/project/badge/Grade/0f4ae5dc6e1140ad855a3d6325d44b35)](https://app.codacy.com/gh/chaotic-society/theoretica/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=chaotic-society/theoretica&amp;utm_campaign=Badge_Grade) [![C++20](https://img.shields.io/badge/C%2B%2B-20-blue.svg)](https://isocpp.org/)
[![PRs Welcome](https://img.shields.io/badge/PRs-welcome-orange.svg)](https://github.com/chaotic-society/theoretica/fork)

> **A C++ math library for scientific computing with a simple and elegant interface.**

Written in modern C++, Theoretica provides a comprenhensive suite of numerical methods designed for **high-performance scientific computing**, without the complexity and steep learning curve typically associated with numerical and performance-critical software.

From simulating complex physical systems, to building machine learning models from scratch, Theoretica offers the speed of C++ with the readability of Python.

If you'd like to join us, to learn or to bring your expertise, make sure to read the [Onboarding Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/ONBOARDING.md) and the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md). Your help is valuable!


## ✨ Why Theoretica?

* **🧠 Elegant & Intuitive API**: Mathematical expressions should look like math, not like machine code.
* **⚡ Blazing Fast**: Built with zero-overhead abstractions to squeeze out every drop of performance from your CPU.
* **🔋 Batteries Included**: Out-of-the-box support for advanced numerical methods, without any dependencies.
* **🖥️ HPC Support**: Growing support for hardware accelleration, making your simulations even faster.
* **🎚️ Embedded Friendly**: Runs anywhere, even with low resources, making it a great choice for embedded systems.

Theoretica is constantly developed and improved with new ideas!


## 🚀 Capabilities
[![Test on Linux](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-linux.yml) [![Test on Windows](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-windows.yml) [![Test on MacOS](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml/badge.svg)](https://github.com/chaotic-society/theoretica/actions/workflows/test-macos.yml)

Theoretica handles the heavy lifting across multiple scientific domains:
* **Simulation & Physics**: Differential equation solvers, Monte Carlo methods, integral quadratures.
* **Numerical Linear Algebra**: Standalone linear algebra algorithms for dense matrices.
* **Automatic Differentiation**: Evaluate exact derivatives without manual analytical calculations or numerical approximations.
* **Optimization & Root Finding**: Built-in routines for finding function extrema and roots, in arbitrary dimensions.
* **Statistics & Error Propagation**: Methods for regression, fits, histograms, data sampling, and built-in error propagation mechanics.
* **Signal Processing**: Standalone Fast Fourier Transform support.
* **File I/O**: First-class support for storing large-scale data sets using CSV and HDF5 file formats, as well as Data Frames.

On each commit, tests are run and documentation is built and deployed online. This ensures that the library works correctly and the documentation is always up-to-date.


## 🦋 Showcase: Simulating Chaos
The following code solves a differential equation, such as the [Lorenz attractor](https://en.wikipedia.org/wiki/Lorenz_system):

```cpp
vec3 f(real t, vec3 v) {

    const real a = 13, b = 20, c = 8/3.0;
    const real x = v[0], y = v[1], z = v[2];

	return {
		a * (y - x),
		x * (b - z) - y,
		x * y - c * z
	};
}

int main() {

    // Solve the system using Runge-Kutta's method
    vec3 v0 = {0.1, 0, 0};
    auto solution = ode::solve_rk4(f, v, 0.0, 50.0);

    // Write the solution directly to file
    std::ofstream file ("attractor.csv");
    file << solution;
}

```

You can find many more examples in the [examples](https://github.com/chaotic-society/theoretica/tree/master/examples) folder.

## 🛠️ Setup

**Theoretica is a header-only library and has no dependencies, so you can include it in your projects straight-away!** To use the library, you can include single headers or use `theoretica.h` which includes all modules, or alternatively include `theoretica_mini.h` which includes only base modules.

You can also compile tests and example programs using Make (`make test` and `make examples`) or using CMake:

```bash
cd build
cmake ..
# Use your chosen compiler
```

## 📖 Documentation
[![Documentation](https://img.shields.io/badge/Doxygen-docs-blue?style=flat&cacheSeconds=https%3A%2F%2Fchaotic-society.github.io%2Ftheoretica%2F&link=https%3A%2F%2Fchaotic-society.github.io%2Ftheoretica%2F)](https://chaotic-society.github.io/theoretica)

The documentation for the project is available [here](https://chaotic-society.github.io/theoretica). The documentation is written using Doxygen syntax alongside the source code and the online version is automatically updated on each commit. The bibliography used during research for the library is listed in the [Bibliography](https://github.com/chaotic-society/theoretica/blob/master/txt/BIBLIOGRAPHY.md). To learn more about the design choices behind the library, you can read the RFC documents in this [folder](https://github.com/chaotic-society/documents/tree/main/specification/theoretica/rfc).


## 🤝 Join Us!

We believe that the best scientific software is built by communities. Whether you are a mathematician, programmer, physicist or student, there is a place for you here.

**How you can make an impact:**
* 🐛 **Report Bugs**: Found an edge case in a solver? Open an issue!
* 💡 **Suggest Features**: Need a specific mathematical distribution or integrator? Let us know.
* ⌨️ **Write Code**: Check out our `good first issue` labels. We are always looking for new solvers, optimizations and features.
* 📖 **Improve Docs**: Help us make Theoretica accessible to scientists and students worldwide.

Have a look at the [Contributing Guide](https://github.com/chaotic-society/theoretica/blob/master/txt/CONTRIBUTING.md) to learn more!

## 📑 License
[![License](https://img.shields.io/github/license/chaotic-society/theoretica)](https://choosealicense.com/licenses/lgpl-3.0/)

The project is currently under the [GNU Lesser General Public License 3.0](https://github.com/chaotic-society/theoretica/blob/master/LICENSE). You may learn more about it [here](https://choosealicense.com/licenses/lgpl-3.0/).

---
