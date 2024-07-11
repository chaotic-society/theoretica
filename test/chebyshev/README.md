# Chebyshev Test
[![Build Status](https://github.com/chaotic-society/chebyshev/actions/workflows/build.yml/badge.svg)](https://github.com/chaotic-society/chebyshev/actions/workflows/build.yml)

Chebyshev is a header-only C++ testing framework designed for testing scientific software and scientific computing libraries. It is part of the larger Theoretica project, a mathematical library that is thoroughly tested using Chebyshev. The framework is composed of three modules: a `prec` module for precision testing, a `benchmark` module for benchmarks, and the `err` module for error checking. Chebyshev provides a robust and flexible way to ensure the accuracy, performance, and reliability of scientific computing applications.

## Features
- **Header-only:** Chebyshev is a header-only library, making it easy to integrate into existing projects without requiring additional dependencies or build steps.

- **Error integral approximation:** Chebyshev provides estimators for the precision of functions over their domain.

- **Modular design:** The framework is composed of three independent modules, allowing developers to use only the components they need.

- **Flexible testing:** Chebyshev provides a range of testing functions and is implemented using templates, making it easy to write custom tests for specific use cases.

- **Many output formats:** The framework supports multiple output formats, such as Markdown and LaTeX tables, and is fully extensible with custom formats.

- **Multi-platform support:** Chebyshev is designed to work on various platforms, including Windows, Linux, and MacOS, and is fully platform-independent.

## Interface
The different modules are contained in their respective namespaces and are initialized through the `<module>::setup()` function and are terminated with the `<module>::terminate()` function, which outputs the results. The behavior of a module may be customized and extended by modifying the fields of the `<module>::state` structure after the module has been initialized. The results of testing are also output to file in CSV or other formats for easy analysis, manipulation and visualization.


### Precision testing
The precision testing module, implemented in the `prec` namespace, is designed to verify the accuracy of scientific computing algorithms. It provides a set of functions to compare the results of different implementations, ensuring that they produce identical or equivalent results within a specified tolerance. In addition to checking single equivalences, Chebyshev implements precision estimation techniques, which consist in estimating error integrals of functions over a certain domain. This generally consists in estimating, either with deterministic quadrature methods or Monte Carlo methods, the following integrals:

$$\epsilon_{mean} = \frac{1}{\mu(\Omega)} \int_\Omega |f(x) - f'(x)| dx$$

$$\epsilon_{rms} = \frac{1}{\mu(\Omega)} \sqrt{\int_\Omega |f(x) - f'(x)|^2 dx}$$

$$\epsilon_{max} = \max_{\Omega} |f(x) - f'(x)|$$

$$\epsilon_{rel} = \frac{\int_\Omega |f(x) - f'(x)| dx}{\int_\Omega |f(x)|dx}$$

The implementation is generalized using templates, making it possible to test quite generic types of functions, from real functions to functions of matrices and vectors or complex numbers. The estimates are computed and the single test cases are validated through a _fail function_, which determines whether the test failed, depending on its results.


### Benchmarks
The `benchmark` module is used to measure the performance of algorithms and functions in general. It provides a set of macros and functions to time and profile the execution of code, allowing developers to optimize their implementations for speed and efficiency. The `benchmark::benchmark()` function works by running the function under consideration for multiple _runs_ and _iterations_, where runs use the same input, while different iterations use different inputs. The average runtime is then computed and registered. The input to feed the function can be fully customized using, for example, randomized input over the domain of the function.


### Error checking
The `err` module makes it possible to test that functions correctly set `errno` or throw exceptions. This is achieved for example by calling the functions with values outside of their domain, checking that they report the error correctly. The functions `err::check_errno` and `err::check_exception()` are used for these type of checks.


### Output customization
The additional `output` module, not directly used for testing, makes it possible to customize the output of the tests, such as which fields to print and how to display them. Customization options are available through the `output::state` structure and are automatically applied to the other modules.


### Randomized tests
The `random` module works in conjunction with the three testing modules to randomize test inputs and provide distribution sampling capabilities for your test units.


## Getting Started
To use Chebyshev, simply include the relevant header file for the module that you need in your project and start writing tests. You can alternatively include the `chebyshev.h` header file which automatically includes all functionalities. The framework is designed to be easy to use, with a minimal learning curve and customization options to shape it according to your needs. This example code sets up precision testing for the "example" test unit and estimates the error over a fictitious function with respect to an exact function:

```c
prec::setup("example", argc, argv);

	// Estimate errors on f on [0, 100]
	prec::estimate("f", f, g, prec::interval(0, 100));

	// Check that two values are equal up to a tolerance
	prec::equals("f", f(1), 1, 0.2);

prec::terminate();
```

The output of this simple code, using default options, is:

```
Starting precision testing of the example module ...

 ┌──────────────────────────────────────────────────────────────────────────────┐
 │         Function │    Mean Err. │     RMS Err. │     Max Err. │       Result │
 ├──────────────────────────────────────────────────────────────────────────────┤
 │                f │      3.3e-12 │      3.5e-12 │      5.1e-12 │         PASS │
 └──────────────────────────────────────────────────────────────────────────────┘

Results have been saved in: example_results

 ┌───────────────────────────────────────────────────────────────┐
 │         Function │   Difference │    Tolerance │       Result │
 ├───────────────────────────────────────────────────────────────┤
 │                f │     -0.0e+00 │      2.0e-01 │         PASS │
 └───────────────────────────────────────────────────────────────┘

Results have been saved in: example_results
Finished testing example
2 total tests, 0 failed (0%)
```
## Setup and Usage
Chebyshev is a header-only library, so there is no need to build or install it separately. Simply include the relevant header files in your project and start using the framework straightaway. Only a compiler with C++14 support is needed to use the framework.


## Contributing
Chebyshev is a collaborative open-source project, and contributions are welcome. If you'd like to contribute to the framework, please submit a pull request.
