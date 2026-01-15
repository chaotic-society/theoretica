# Testing

The [Chebyshev](https://github.com/chaotic-society/chebyshev) testing framework is used to validate Theoretica's numerical methods. In particular, the `prec` module is used to estimate the precision of the methods, while the `benchmark` module is used to benchmark performance-critical code.

## Directory Structure
All test programs are stored inside the `test` directory and Make may be used to execute them, using `make test` for precision testing and `make benchmark` for benchmarking. Precision testing programs are stored in `test/prec`, with implementation files organized by module and called `test_modulename.cpp`, while benchmarks are inside the `test/benchmark` folder and are called `benchmark_modulename.cpp`.

### Precision Testing
The `prec` module provides two main functions, `estimate()` and `equals()`, implemented inside the `prec_context` class, which handles all internal state. The `estimate()` function uses an estimator function to estimate error integrals over a domain of a certain function with respect to an exact or better implementation (or even with respect to a certain property of the function). Chebyshev provides some built-in estimators for common applications, such as real functions using quadrature of integrals or Monte Carlo methods. On the other hand, `equals()` compares two values, checking that their distance is smaller than the tolerance (with respect to an arbitrary distance function). The context is created by initializing a `prec_context` class or using `prec::make_context()`:

```cpp
auto ctx = prec::make_context("moduleName");

	// Test cases ...

// Optional
ctx.terminate();
```

The settings for the module are available by modifying the `ctx.settings` structure directly. For example, the following code sets the output files:

```cpp
ctx.settings.outputFiles = { "output.csv" };
```

By default, the output format is CSV, as is generally used for Theoretica's test output, but it may be changed, and custom formats are supported. The `estimate_options` and `equation_options` structures are used to store additional options which may be used between multiple test cases. This code shows how to use `estimate()` with a custom estimator:

```cpp
auto options = prec::estimate_options<real, real>(
	prec::interval(0.0, 1.0),			// Interval of estimation
	prec::estimator::quadrature1D(),	// 1D integral quadrature
	1E-04								// Tolerance
);

ctx.estimate(
	"myFunctionName", approxFunction, exactFunction, options
);
```

The template parameters to `estimate_options` are the return and argument types of the function under test, while the constructor arguments are the name of the test case (usually the function name), the function under test and a reference function, either exact or more accurate (sometimes it's easier to test over known equations with expected function equal to zero or a constant). The last argument provides the options for the test case, with the first argument being an interval or list of intervals which determine the domain, the second being the estimator to use and the last the tolerance value for the tests. 

The following code uses the `equals()` function to compute the distance between two values, in this case two polynomials, using a custom provided distance function:

```cpp
auto options = prec::equation_options<polynomial<>>(
	1E-08, distance_polyn
);

ctx.equals(
	"polynomial", P1, P2, options
);
```

The template parameter for `equation_options` is the type of variable which will be compared, while `distance_polyn` is a custom distance function taking in two polynomials and computing their distance over the coefficients (the first argument to the constructor is the tolerance value). If a distance function is not provided, it defaults to NaN for custom types (which makes the test fail) and `abs()` for floating point numbers.
