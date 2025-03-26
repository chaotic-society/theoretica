# Testing

The [Chebyshev](https://github.com/chaotic-society/chebyshev) testing framework is used to validate Theoretica's numerical methods. In particular, the `prec` module is used to estimate the precision of the methods, while the `benchmark` module is used to benchmark performance-critical code.

## Directory Structure
All test programs are stored inside the `test` directory and Make may be used to execute them, using `make test` for precision testing and `make benchmark` for benchmarking. Precision testing programs are stored in `test/prec`, with implementation files organized by module and called `test_modulename.cpp`, while benchmarks are inside the `test/benchmark` folder and are called `benchmark_modulename.cpp`.

### Precision Testing
The `prec` module provides two main functions, `prec::estimate` and `prec::equals`. The `prec::estimate` function uses an estimator function to estimate error integrals (or sums) over a domain of a certain function with respect to an exact or better implementation. Chebyshev provides some built-in estimators for common applications, such as real functions using quadrature of integrals or Monte Carlo methods. On the other hand, `prec::equals` compares two values, checking that their distance is smaller than the tolerance (with respect to an arbitrary distance function). The module is initialized using `prec::setup("module name")` and terminated using `prec::terminate()`:

```cpp
prec::setup("module name");

	//...

prec::terminate();
```

The test cases are written between the two directives. The settings for the module are available by modifying the `prec::settings` structure directly. For example, the following code sets an output file:

```cpp
prec::settings.outputFiles = { "output.csv" };
```

By default, the output format is CSV, as is generally used for Theoretica's test output. The `estimate_options` and `equation_options` structures are used to store additional options which may be used between many different test cases. This code shows how to use `prec::estimate` with a custom estimator:

```cpp
auto options = prec::estimate_options<real, real>(
	prec::interval(0.0, 1.0),
	prec::estimator::quadrature1D(),
	1E-04
);

prec::estimate(
	"derivative",
	approxFunction, exactFunction, options
);
```

The template parameters to `estimate_options` are the return and argument types of the function under test, while the constructor arguments are the name of the test case (usually the function name), the function under test and a reference function, either exact or more accurate (sometimes it's easier to test over known equations with expected function equal to zero or a constant). The last argument provides the options for the test case, with the first argument being an interval or list of intervals which determine the domain, the second being the estimator to use and the last the tolerance value for the tests. The following code uses the `prec::equals` function to compute the distance between two values, in this case two polynomials, using a custom provided distance function:

```cpp
auto options = prec::equation_options<polynomial<real>>(
	1E-08, distance_polyn
);

prec::equals(
	"polynomial",
	P1, P2, options
);
```

The template parameter for `equation_options` is the type of variable which will be compared, while `distance_polyn` is a custom distance function taking in two polynomials and computing their distance over the coefficients (the first argument to the constructor is the tolerance value). If a distance function is not provided, the ordering `operator<` is used to compute the absolute value of any ordered type.
