# 🧪 Testing Guide

Testing is the backbone of any reliable scientific computing library.  When dealing with numerical methods even small errors can propagate and lead to significant inaccuracies. That's why we take testing seriously, and why your contributions to our test suite are so valuable.

Unlike typical software where a function either works or doesn't, numerical methods can have a spectrum of outcomes. A square root function might return `1.9999` instead of `2.0`, and that's often acceptable, but how do we know when the error is too large? This is where precision testing comes in, with features such as:

- **Estimating numerical errors** across entire domains of input values
- **Comparing approximate results** against known exact solutions
- **Measuring different error metrics** (mean error, maximum error, RMS error)
- **Benchmarking performance** to ensure implementations are efficient

## Chebyshev Framework

We use [Chebyshev](https://github.com/chaotic-society/chebyshev), our custom-built testing framework, to validate Theoretica's numerical methods. It provides three main modules:

| Module | Purpose | When to Use |
|--------|---------|-------------|
| **`prec`** | Precision testing | Verifying numerical accuracy of mathematical functions |
| **`benchmark`** | Performance measurement | Measuring execution time |
| **`err`** | Error handling | Testing error conditions and assertions |

For most contributions, you'll primarily work with the **`prec`** module.

## Directory Structure

All test programs live inside the `test` directory. Here's how everything is organized:

```
test/
├── chebyshev/              # The Chebyshev framework (git submodule)
├── prec/                   # Precision testing programs
│   ├── test_*.cpp			# Test unit
│   └── test_template       # Template for new test units
└── benchmark/              # Performance benchmarks
	├── benchmark_*.cpp		# Benchmark unit
	└── benchmark_template	# Template for new benchmarks
```

Each test file in `test/prec/` corresponds to a module in the library. For example, `test_calculus.cpp` tests functions from the `calculus` module (derivatives, integrals, ODE solvers), while `test_complex.cpp` tests complex number arithmetic.

## Running Tests

You can use Make to compile and run all tests:

```bash
make all		# Build everything and run all tests
make test		# Run only precision tests
make benchmark	# Run only benchmarks
```

Or equivalently CMake:

```bash
cd build
cmake ..
make test
make benchmark
```

If you create a new test unit, these commands will automatically run your test too! To run your test unit using the build system, you can use `make test_modulename`.

## Precision Testing

The `prec` module provides two essential functions: `ctx.estimate()` for testing functions over continuous domains, and `ctx.equals()` for comparing specific values. 

### Creating a Test Unit

Every test file begins by creating a **context**, an object that manages test execution, collects results, and handles output. The specific type for accuracy testing is `prec_context`. Think of it as the "test runner" that orchestrates everything: 

```cpp
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;

int main(int argc, char const *argv[]) {
    
    auto ctx = prec::make_context("moduleName", argc, argv);
    ctx.settings.outputFiles = { "test/prec/prec_moduleName.csv" };

    // Your test cases go here...

    // The context automatically cleans up when it goes out of scope,
    // but you can explicitly terminate if needed: 
    ctx.terminate();
}
```

The context object `ctx` provides methods for running tests and configuring behavior. You can customize various settings: 

```cpp
// Configure which columns appear in the output
ctx.settings.estimateColumns = {
    "name", "meanErr", "rmsErr", "maxErr", "tolerance", "failed"
};

// Enable or disable multithreading
ctx.settings.multithreading = false;
```

### The `estimate()` Function

The `ctx.estimate()` function is used for testing continuous mathematical functions. Instead of checking just one or two values, it evaluates your function across an entire domain and estimates error integrals and other heuristics.

**When to use `estimate()`:**
- Testing continuous mathematical functions
- Checking that an approximation is accurate across its intended domain

**How it works:**
The function takes your approximate implementation and an exact (reference) function, then compares them at many points within a specified interval. It computes several error metrics, such as mean error or maximum error.

The test domain is parametrized as a product of sub-intervals of the real numbers $\Omega \subset \mathbb{R}^n$. For example, a function of two variables `x` and `y` may be tested over a domain $[0, 1] \times [3, 5]$, which in code is expressed as `{ prec::interval(0, 1), prec::interval(3, 5) }`. The specific interpretation of the domain is left to the chosen **estimator**.

Here's a basic example:

```cpp
// Custom options for our test case
auto options = prec::estimate_options<real, real>(
    prec::interval(0.0, 2 * PI),		// Domain: [0, 2 pi]
    prec::estimator::quadrature1D(),	// Estimation method
    1E-10								// Tolerance
);

ctx.estimate(
	"my_sin",			// Test case name
	my_sin,				// Our implementation
	std::sin,           // Reference implementation
	options				// Custom options
);
```

Error metrics are estimated using an object called **estimator**, which takes in the approximate function, the exact function and other options such as the domain, and uses a certain numerical method. For example, `prec::estimator::quadrature1D()` estimates error integrals for real functions in 1 variable using simple quadrature methods, while `prec::estimator::montecarlo(n)` estimates error integrals for multivariate functions using Monte Carlo methods. Depending on the functions being tested, it may be necessary to write custom estimators. In this case, don't be afraid to reach out to maintainers for help!

**Real-world example** testing numerical differentiation:

```cpp

real f(real x) {
	// ...
}

real Df(real x) {
	// ...
}

auto deriv_opt = prec::estimate_options<real, real>(
    prec::interval(0, 1),				// Test for x in [0, 1]
    prec::estimator::quadrature1D(),	// Use 1D quadrature
    1E-04,								// Tolerance is 10^-4
	1000								// Use 1000 points
);

ctx.estimate("deriv_central",
    [](real x) { return deriv_central(f, x, 10E-8); },
    Df, deriv_opt
);
```

The `estimate_options` template takes two or more type parameters:
- The **return type** of the function being tested
- The **argument types** of the function being tested

For a function `int f(real x)`, you'd use `estimate_options<int, real>`, while for a function `int g(real x, real y)` you would use `estimate_options<int, real, real>`.

### The `equals()` Function
While `estimate()` tests over ranges, `equals()` compares two specific values. This is ideal for discrete tests or cases where you have known exact results.

**When to use `equals()`:**
- Verifying edge cases
- Checking specific known values (e.g., `sqrt(4) = 2`)

```cpp
// Simple equality check
ctx.equals("sqrt(4) == 2", th::sqrt(4.0), 2.0);
```

For types more complex than numbers (like polynomials, matrices, or vectors) you may need to provide a custom distance function that measures how "far apart" two values are: 

```cpp
// Custom distance function for polynomials
prec_t distance_polyn(const polynomial<real>& p1, const polynomial<real>& p2) {
    const polynomial<real> d = p1 - p2;
    real r = -inf();

    for (size_t i = 0; i < d.size(); ++i)
        r = max(r, std::abs(d[i]));

    return r;
}

// Use the custom distance in options
auto opt = prec::equation_options<polynomial<real>>(
    1E-08,          // Tolerance
    distance_polyn  // Custom distance function
);

polynomial<real> p = {1. 0, 1.0, 1.0};		// 1 + x + x²
polynomial<real> expected = {1.0, 2.0};		// 1 + 2x (the derivative)

ctx.equals("deriv(P)", deriv(p), expected, opt);
```

A distance function takes in two values of a given type and results a non-negative `prec_t` number representing their distance.

### Choosing Between `estimate()` and `equals()`

| Scenario | Use `estimate()` | Use `equals()` |
|----------|------------------|----------------|
| Testing `sin(x)` over [0, 2π] | ✅ | |
| Checking that `sqrt(4) = 2` | | ✅ |
| Verifying numerical integration | ✅ | |
| Comparing two polynomials | | ✅ |
| Testing iterative solvers | | ✅ |
| Measuring derivative approximation error | ✅ | |

## Benchmarking
Benchmarks measure how fast a certain implementation runs. They live in `test/benchmark/` and use Chebyshev's `benchmark` module. Here's a simple example:

```cpp
auto ctx = benchmark::make_context("my_module", argc, argv);
ctx.settings.outputFiles = { "test/benchmark/benchmark_my_module.csv" };

auto opt = benchmark::benchmark_options<real>(
	10,			// Number of runs (same inputs)
	100000,		// Number of iterations per run
    benchmark::generator::uniform1D(0, +1E+06) // Random number generator
);

ctx.benchmark("my_function", my_function, opt);
```

**Benchmarking Tips:**
- Disable multithreading when benchmarking PRNGs or other stateful operations
- Use enough iterations to get stable measurements (at least 10k for fast functions)
- Test with realistic input ranges that match actual use cases
- Compare against baseline implementations when optimizing

## Best Practices

### Naming Conventions
- Precision test files:  `test_<module>.cpp`
- Benchmark files: `benchmark_<module>.cpp`
- Output CSV files should match the source file name (except for extension)

### Choosing Tolerances
Tolerance values depend on the function you're testing, but generally speaking, Theoretica aims for a baseline 1E-08 max error. Depending on the numerical method, this may be too small (the method has limited convergence, such as `th::deriv_forward()`), or could be made smaller (the function is very accurate, such as `th::sqrt()`).

### Writing Good Tests

1. **Test bulk domain**: Ensure functions work correctly under normal conditions
2. **Test edge cases**: Zero, negative values, very large numbers, infinities
3. **Test numerical stability**: Operations near the limits of floating-point precision
4. **Document your tests**: Add comments explaining what each test verifies, if unclear

## Getting Help

- **Look at existing tests** in `test/prec/` for examples
- **Check [Chebyshev documentation](https://github.com/chaotic-society/chebyshev)**
- **Ask a maintainer**: we're happy to help!

Remember, imperfect tests are better than no tests. Start simple, and we can refine together.

### **Thank you for helping us build a more robust library!**
