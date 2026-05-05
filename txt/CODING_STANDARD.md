# Theoretica Coding Standard

Welcome! If you're contributing to Theoretica, this guide will help you write code that's consistent with the rest of the library. Following these standards improves code quality, readability, and makes collaboration smoother for everyone.

**Thank you for your contribution!** 🎉

---

## Table of Contents

1. [Project Philosophy](#project-philosophy)
2. [Namespace & Scope](#namespace--scope)
3. [Naming Conventions](#naming-conventions)
4. [Code Formatting](#code-formatting)
5. [Documentation](#documentation)
6. [Header Files & Inline Functions](#header-files--inline-functions)
7. [Templates](#templates)
8. [Error Handling](#error-handling)
9. [Testing](#testing)
10. [Performance Considerations](#performance-considerations)
11. [Examples of Good Code](#examples-of-good-code)
12. [Common Pitfalls](#common-pitfalls)
13. [Before Submitting](#before-submitting)

---

## Project Philosophy

Theoretica aims to be:
- **Mathematically elegant:** Code should mirror mathematical notation where possible
- **Simple and readable:** Favor clarity over cleverness
- **Performant:** Scientific computing demands efficiency
- **Functional-style:** Prefer pure functions without side effects when possible
- **Type-safe:** Leverage C++ type system for correctness

Keep these principles in mind when writing code!

---

## Namespace & Scope

### Namespaces

**All functions and classes MUST be declared inside the `theoretica` namespace:**

```cpp
namespace theoretica {

	// Your code here

}
```

**Sub-namespaces** should only be used when grouping related functionality. Examples:
- `theoretica::algebra` for algebraic operations
- `theoretica::taylor` for Taylor series expansions
- `theoretica::stats` for statistical functions

**Example:**
```cpp
namespace theoretica {
namespace taylor {
	
	inline polynomial<real> expand_linear(/* ... */) {
		// Implementation
	}
		
}
}
```

### Global Variables & Constants

- **Global variables:** ❌ Avoid unless absolutely necessary
- **Global constants:** ✅ Acceptable and encouraged

**Global constants must use SCREAMING_SNAKE_CASE:**

```cpp
namespace theoretica {
		
	/// Gravitational constant in SI units
	constexpr real GRAVITATIONAL_CONSTANT = 6.67430e-11;

}
```

**Configuration constants** follow the pattern `THEORETICA_<CONTEXT>_<NAME>`:

```cpp
#define THEORETICA_NO_PRINT		// Disables printing functionality
#define THEORETICA_VECTOR_SIZE	// Maximum static vector size
```

---

## Naming Conventions

### Functions

**Use snake_case for all function names:**

```cpp
inline real square_root(real x);
inline mat<real> decompose_lu(const mat<real>& A);
```

**Naming pattern:** `action_method_modifier`

| Component | Description | Examples |
|-----------|-------------|----------|
| **action** | What the function does | `compute`, `find`, `solve`, `decompose` |
| **method** | What algorithm/method is used | `lu`, `qr`, `newton`, `bisection` |
| **modifier** | Additional specifics (optional) | `inplace`, `approx`, `fast` |

**Examples:**
```cpp
// Good naming
inline void decompose_lu_inplace(mat<real>& A);
inline real integral_simpson(std::function<real(real)> f, real a, real b);

// Avoid
inline void doLU(mat<real>& A);					// Wrong case
inline real super_fast_integrator(/* ... */);	// Not following pattern
```

### Classes

**Use snake_case for class names:**

```cpp
class complex;
class dual;
class mat_iterator;
```

**Class names should:**
- Represent the mathematical concept they encapsulate
- Be concise but descriptive
- Avoid abbreviations unless universally understood

**Examples:**
```cpp
// Mathematical structures
class complex { /* ... */ };
class quaternion { /* ... */ };
class dual { /* ... */ };
class polynomial { /* ... */ };

// Algorithms/Utilities
class mat_iterator { /* ... */ };
class PRNG { /* ... */ };	// Acceptable: well-known acronym
```

### Constants

**Use SCREAMING_SNAKE_CASE:**

```cpp
constexpr real PI = 3.14159265358979323846;
constexpr real EULER_GAMMA = 0.57721566490153286060;
constexpr real MACH_EPSILON = 1e-16;
```

### Macros

**Use macros sparingly.** When necessary, follow `TH_<NAME>`:

```cpp
#define TH_MATH_ERROR(func, arg, type) \
		// Error handling macro

#define TH_NO_THROW noexcept
```

**Prefer `constexpr` or `inline` functions over macros when possible.**

### Template Parameters

**Use CamelCase for template type parameters:**

```cpp
template<typename Type>
class vector { /* ... */ };

template<typename ReturnType, typename InputType>
inline ReturnType convert(InputType x);

template<typename Matrix, typename ElementType>
class mat_iterator { /* ... */ };
```

**Examples from the codebase:**
```cpp
// From mat.h
template<typename Matrix, typename ReturnType = matrix_element_t<Matrix>&>
class mat_iterator { /* ... */ };

// From multidual.h
template<unsigned int N = 0>
class multidual { /* ... */ };
```

### Member Variables

**Use simple lowercase names for class members:**

```cpp
class complex {
public:
	real a;	// Real part
	real b;	// Imaginary part
};
```

**For iterators and more complex classes, use descriptive snake_case:**

```cpp
class mat_iterator {
private:
	Matrix& matrix;
	size_t row;
	size_t col;
};
```

---

## Code Formatting

### Indentation

**Use TABS for indentation (not spaces):**

```cpp
namespace theoretica {
	
	inline real square(real x) {
		return x * x;
	}
	
	class complex {
	public:
		
		complex(real r, real i) : a(r), b(i) {}
		
		inline real norm() const {
			return sqrt(a * a + b * b);
		}
	};
	
}
```

### Spacing

**Add a space after commas:**

```cpp
// Good
function(a, b, c);
template<typename T, typename U>
std::array<real, 3> values = {1, 2, 3};

// Bad
function(a,b,c);
template<typename T,typename U>
```

**Spacing around operators (recommended but flexible):**

```cpp
// Preferred
x = a + b * c;
y = (a + b) / (c - d);

// Also acceptable
x=a+b*c;	// If it improves readability in math-heavy contexts
```

### Braces

**Place opening braces on the same line:**

```cpp
// Good
if (condition) {
	// code
}

// Also acceptable for single-line functions
inline real square(real x) { return x * x; }
```

### Line Length

**No strict limit, but aim for readability:**
- Prefer line breaks for complex expressions
- Break long function signatures logically

```cpp
// Good
inline polynomial<real> expand_quadratic(
	Dual2Function f,
	real x0 = 0
) {
	// Implementation
}

// Also acceptable if not too long
inline polynomial<real> expand_linear(DualFunction f, real x0 = 0) {
	// Implementation
}
```

---

## Documentation

### Doxygen Format

**All public functions and classes MUST be documented using Doxygen:**

**Basic function documentation:**

```cpp
/// Compute the square root of a number
///
/// @param x A non-negative real number
/// @return The square root of x
///
/// Throws an error if x is negative.
inline real sqrt(real x) {
	// Implementation
}
```

**Detailed function documentation:**

```cpp
/// Compute the LU decomposition of a matrix in-place
///
/// @param A The matrix to decompose (modified in-place)
/// @param L Output lower triangular matrix (optional)
/// @param U Output upper triangular matrix (optional)
/// @return True if decomposition succeeded, false otherwise
///
/// This function uses Doolittle's algorithm to compute the LU decomposition
/// such that \f$A = LU\f$ where L is lower triangular and U is upper triangular.
///
/// @note The input matrix A is modified in-place to store the result
/// @warning This function requires A to be square and non-singular
///
/// **Error conditions:**
/// - If A is singular, the function returns false
/// - If A is not square, behavior is undefined
///
/// **Example:**
/// ```cpp
/// mat<real> A = {{4, 3}, {6, 3}};
/// mat<real> L, U;
/// decompose_lu_inplace(A, L, U);
/// // A is now modified; L and U contain the factors
/// ```
inline bool decompose_lu_inplace(
		mat<real>& A,
		mat<real>* L = nullptr,
		mat<real>* U = nullptr
) {
		// Implementation
}
```

**Class documentation:**

```cpp
/// @class complex
/// Complex number in algebraic form \f$a + ib\f$
///
/// This class represents complex numbers and provides
/// arithmetic operations, trigonometric functions, and
/// conversions between different representations.
///
/// @note Template parameter Type defaults to `real`
template<typename Type = real>
class complex {
	// Implementation
};
```

### Mathematical Formulas

**Use LaTeX notation between `\f$` and `\f$`:**

```cpp
/// Compute the Euclidean norm
///
/// The norm is defined as \f$\|x\| = \sqrt{x_1^2 + x_2^2 + \ldots + x_n^2}\f$
///
/// @param v A vector
/// @return The Euclidean norm of v
inline real norm(const vec<real>& v);
```

**For display-style formulas, use `\f[` and `\f]`:**

```cpp
/// Compute the definite integral using Simpson's rule
///
/// The approximation is given by:
/// \f[
///	 \int_a^b f(x) dx \approx \frac{h}{3}\left[f(a) + 4f\left(\frac{a+b}{2}\right) + f(b)\right]
/// \f]
/// where \f$h = b - a\f$.
inline real integral_simpson(std::function<real(real)> f, real a, real b);
```

### Documentation Tags

Common Doxygen tags to use:

| Tag | Purpose | Example |
|-----|---------|---------|
| `@param` | Describe a parameter | `@param x The input value` |
| `@tparam` | Describe a template argument | `@tparam Type The input type` |
| `@return` | Describe return value | `@return The computed result` |
| `@note` | Additional information | `@note This modifies the input` |
| `@warning` | Important warnings | `@warning x must be positive` |
| `@see` | Reference related functions | `@see sqrt, square` |
| `@todo` | Mark incomplete features | `@todo Add error handling` |

---

## Header Files & Inline Functions

### Inline Declarations

**All functions except class constructors SHOULD be `inline` when defined in headers:**

```cpp
// Good
inline real square(real x) {
	return x * x;
}

// Good - constructor doesn't need inline
class complex {
public:
	complex(real r, real i) : a(r), b(i) {}
	
	// Member functions should be inline
	inline real norm() const {
		return sqrt(a * a + b * b);
	}
};
```

**Why?** Header-only libraries require `inline` to avoid multiple definition errors.

### Include Guards

**Use `#ifndef` guards (not `#pragma once`):**

```cpp
#ifndef THEORETICA_FILENAME_H
#define THEORETICA_FILENAME_H

// Header contents

#endif
```

**Naming pattern:** `THEORETICA_<FILENAME>_H`, or optionally `THEORETICA_<MODULE>_<FILENAME>_H` when clashing names are possible.

**Examples:**
```cpp
#ifndef THEORETICA_COMPLEX_H
#ifndef THEORETICA_MATRIX_H
#ifndef THEORETICA_INTERP_POLYNOMIAL_H
```

### Using Directives in Headers

**❌ NEVER use `using namespace` in header files:**

```cpp
// BAD - pollutes global namespace for all users
using namespace std;
using namespace theoretica;
```

**✅ It's OK in implementation files (.cpp), test files, and examples:**

```cpp
// examples/random_walk.cpp
#include "theoretica.h"
using namespace th;	// OK in .cpp files

int main() {
	// ...
}
```

### Conditional Compilation

**For optional features, use guards:**

```cpp
#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

// Later in code:
#ifndef THEORETICA_NO_PRINT
	inline std::string to_string() const {
		// Implementation
	}
#endif
```

---

## Templates

### Template Syntax

**Template parameters should use CamelCase:**

```cpp
template<typename Type>
class vector { /* ... */ };

template<typename InputType, typename OutputType>
inline OutputType convert(InputType x);
```

### Template Constraints (C++14)

**Use `std::enable_if` or SFINAE for template constraints when needed:**

```cpp
// Type traits example
template<typename Matrix>
using matrix_element_t = /* ... */;

// SFINAE example (if needed)
template<typename T>
typename std::enable_if<std::is_floating_point<T>::value, T>::type
safe_sqrt(T x) {
		return sqrt(x);
}
```

The header `core_traits.h` contains useful type traits for Theoretica.

**Note:** In the future (C++20), we may partially migrate to concepts.

### Default Template Parameters

**Prefer defaults that match common use cases:**

```cpp
template<typename Type = real>
class complex { /* ... */ };

template<unsigned int N = 0>	// 0 = dynamic size
class multidual { /* ... */ };
```

---

## Error Handling

### Error Macros

**Use `TH_MATH_ERROR` for mathematical errors:**

```cpp
inline dual inverse() const {
	
	if (a == 0) {
		TH_MATH_ERROR("dual::inverse", 0, MathError::DivByZero);
		return dual(nan(), nan());
	}
	
	return dual(1.0 / a, -b / square(a));
}
```

### Error Types

Common error types in Theoretica:
- `MathError::DivByZero` - Division by zero
- `MathError::InvalidArgument` - Invalid input
- `MathError::OutOfRange` - Index out of bounds
- `MathError::NoConvergence` - Iterative method failed to converge

### Return NaN for Invalid Operations

**For numerical errors, return `nan()` when appropriate:**

```cpp
inline real safe_log(real x) {

	if (x <= 0) {
		TH_MATH_ERROR("safe_log", x, MathError::InvalidArgument);
		return nan();
	}

	return log(x);
}
```

### Assertions vs. Runtime Checks

- **Runtime checks:** For user input that might be invalid
- **Assertions:** For internal consistency checks during development

---

## Testing

### Writing Tests

**Every new feature MUST include tests.** Tests go in `test/prec/test_<module>.cpp`.

**Use the Chebyshev framework:**

```cpp
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;

int main(int argc, char const *argv[]) {
		
	prec_context ctx;
	
	// Test a single value
	ctx.equals("sqrt(4)", sqrt(4.0), 2.0);
	
	// Test over a range
	ctx.estimate(
		"square",
		[](real x) { return square(x); },
		[](real x) { return x * x; },
		prec::interval(0, 10)
	);
	
	return ctx.results();
}
```

### Test File Naming

- **Precision tests:** `test/prec/test_<module>.cpp`
- **Benchmarks:** `test/benchmark/benchmark_<module>.cpp`

### What to Test

**At minimum, test:**
1. **Correctness:** Does it give the right answer?
2. **Edge cases:** What about zero, negative numbers, infinity, NaN?
3. **Boundary conditions:** What about very large or very small inputs?
4. **Known values:** Compare against published tables or exact solutions

**Example test structure:**

```cpp
// Test basic functionality
ctx.equals("complex::operator+()", (complex(1, 2) + complex(3, 4)).Re(), 4.0);

// Test edge case
ctx.equals("sqrt(0)", sqrt(0.0), 0.0);

// Test precision over range
ctx.estimate(
	"sin",
	[](real x) { return th::sin(x); },
	[](real x) { return std::sin(x); },
	prec::interval(0, 2 * PI)
);
```

### Documentation in Tests

**If testing is complex, add comments explaining:**
- What reference values are from (book, paper, calculator)
- Any known limitations
- Why certain tolerances are used

---

## Performance Considerations

### Prefer Pure Functions

**Pure functions** (no side effects, deterministic) should be preferred:

```cpp
// Good - pure function
inline real square(real x) {
	return x * x;
}

// Avoid when possible - side effects
inline void square_inplace(real& x) {
	x = x * x;
}
```

**Exception:** In-place operations are acceptable when performance-critical and clearly marked with `_inplace` suffix.

### Algorithm Complexity

**Document algorithmic complexity for non-trivial functions:**

```cpp
/// Sort a vector using quicksort
///
/// @param v The vector to sort
/// @return A sorted copy of the vector
///
/// **Complexity:** Average O(n log n), worst case O(n²)
inline vec<real> sort_quick(const vec<real>& v);
```

### Avoid Unnecessary Copies

**Use const references for input parameters:**

```cpp
// Good
inline real norm(const vec<real>& v);

// Bad - unnecessary copy
inline real norm(vec<real> v);
```

### Minimize Allocations

**Reuse memory when possible, especially in loops:**

```cpp
// Good
vec<real> result(n);
for (size_t i = 0; i < n; ++i) {
	result[i] = compute(i);
}

// Avoid
vec<real> result;
for (size_t i = 0; i < n; ++i) {
	result.append(compute(i));	// Repeated reallocations
}
```

---

## Examples of Good Code

### Example 1: Simple Function

```cpp
/// Compute the square of a real number
///
/// @param x A real number
/// @return The value of \f$x^2\f$
inline real square(real x) {
	return x * x;
}
```

### Example 2: Class with Methods

```cpp
/// @class complex
/// Complex number in algebraic form \f$a + ib\f$
template<typename Type = real>
class complex {
public:

	/// Real part
	Type a;

	/// Imaginary part
	Type b;


	/// Construct a complex number from real and imaginary parts
	complex(Type real_part, Type imag_part)
		: a(real_part), b(imag_part) {}


	/// Compute the norm of the complex number
	/// @return The value \f$\sqrt{a^2 + b^2}\f$
	inline Type norm() const {
		return sqrt(a * a + b * b);
	}


	/// Compute the complex conjugate
	/// @return The conjugate \f$a - ib\f$
	inline complex conjugate() const {
		return complex(a, -b);
	}
};
```

### Example 3: Template Function

```cpp
/// Compute the Euclidean distance between two vectors
///
/// @param v1 First vector
/// @param v2 Second vector
/// @return The distance \f$\|v_1 - v_2\|\f$
///
/// @warning The vectors must have the same size
template<typename VectorType>
inline real distance(const VectorType& v1, const VectorType& v2) {
	// ...
}
```

---

## Common Pitfalls

### ❌ Don't: Use `using namespace` in Headers

```cpp
// BAD
#ifndef THEORETICA_MODULE_H
#define THEORETICA_MODULE_H

using namespace std;	// ❌ Pollutes global namespace

// ...
#endif
```

### ❌ Don't: Forget `inline` for Functions in Headers

```cpp
// BAD - will cause linker errors
real square(real x) {	// ❌ Missing 'inline'
	return x * x;
}

// GOOD
inline real square(real x) {
	return x * x;
}
```

### ❌ Don't: Use GOTO

```cpp
// BAD
inline real compute(real x) {
	if (x < 0)
		goto error;	// ❌ Never use GOTO
	
	return sqrt(x);
	
error:
	return nan();
}

// GOOD
inline real compute(real x) {
	if (x < 0) {
		TH_MATH_ERROR("compute", x, MathError::InvalidArgument);
		return nan();
	}
	return sqrt(x);
}
```

### ❌ Don't: Ignore Const Correctness

```cpp
// BAD
real norm(vec<real>& v) {	// ❌ Should be const&
	// ...
}

// GOOD
real norm(const vec<real>& v) {
	// ...
}
```

### ❌ Don't: Use Ambiguous Abbreviations

```cpp
// BAD
inline mat decomp_lu(mat A);	// ❌ Abbreviation unclear
inline real calc(real x);		// ❌ Too generic

// GOOD
inline mat decompose_lu(mat A);
inline real sqrt(real x);
```

---

## Before Submitting

### Pre-Submission Checklist

Before you submit a pull request, please ensure:

- [ ] **Code compiles** without warnings (run `make examples`)
- [ ] **Follows naming conventions** (snake_case for functions/classes, CamelCase for templates)
- [ ] **All functions are `inline`** (except constructors)
- [ ] **Documented with Doxygen** (`@param`, `@return`, description)
- [ ] **Inside `theoretica` namespace**
- [ ] **No `using namespace` in headers**
- [ ] **Includes test cases** in `test/prec/test_<module>.cpp`
- [ ] **Tests pass** (run `make test`)
- [ ] **Tabs for indentation** (not spaces)
- [ ] **Spaces after commas**

### Testing Your Code

```bash
# Build and run all tests
make all

# Run only your test
make test_<module>

# Run benchmarks (if applicable)
make benchmark
```

### Asking for Help

Not sure about something? **Don't hesitate to ask!**

- Open a GitHub Discussion
- Comment on your Pull Request draft
- Reference this guide in your questions

We're here to help, especially for first-time contributors!

---

## Quick Reference

### Function Naming Pattern

```
action_method_modifier
   ↓      ↓     ↓
decompose_lu_inplace
```

### Case Styles Quick Guide

| Element | Style | Example |
|---------|-------|---------|
| Functions | `snake_case` | `decompose_lu` |
| Classes | `snake_case` | `complex`, `dual` |
| Templates | `CamelCase` | `Type`, `Matrix` |
| Constants | `SCREAMING_SNAKE_CASE` | `PI`, `EULER_GAMMA` |
| Macros | `TH_<NAME>` | `TH_MATH_ERROR` |
| Members | `lowercase` or `snake_case` | `a`, `b`, `row_count` |

### Documentation Template

```cpp
/// Brief one-line description
///
/// @param param1 Description of first parameter
/// @param param2 Description of second parameter
/// @return Description of return value
///
/// Longer description with details about algorithm,
/// special cases, and mathematical formulas like \f$x^2\f$.
///
/// @note Additional notes
/// @warning Important warnings
/// @see related_function
inline return_type function_name(type1 param1, type2 param2) {
	// ...
}
```

---

## Additional Resources

- **Testing Guide:** See `txt/TESTING.md` for detailed testing information
- **Chebyshev Framework:** https://github.com/chaotic-society/chebyshev
- **Doxygen Documentation:** https://www.doxygen.nl/manual/docblocks.html
- **Project Repository:** https://github.com/chaotic-society/theoretica

---

## Questions?

If anything in this guide is unclear or you have suggestions for improvement, please:

1. Open a GitHub Discussion
2. File an issue with the `documentation` label
3. Submit a PR to improve this document

Happy coding!

---

*This coding standard is a living document and may be updated as the project evolves.*
