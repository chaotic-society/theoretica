# Getting started
Here you can get an overall view of how to use the library for basic use cases, as well as what kinds of mathematical structures are available. The complete documentation is available at [Theoretica Docs](http://chaotic-society.github.io/theoretica).

## Using the library

Theoretica is a **header-only** library, so you only need to include the header files that you need to directly use them. The `theoretica.h` header includes all header files for you (if `THEORETICA_INCLUDE_BASE` is defined, only generic headers such as linear algebra and complex numbers will be included).

```cpp
#include "theoretica.h"
using namespace th;


int main() {
	/// ...
}

```

All functions and classes are defined inside the `theoretica` namespace, with the alias `th`. The `constants.h` header defines common types and constants used throughout the library, such as the `real` type (alias for a floating-point type) and `PI`.

## Linear algebra
The linear algebra module is `algebra`, which defines vector and matrix classes as well as distances, norms and vectorized operations. The vector and matrix classes are divided between *statically* and *dynamically* allocated types. Vectors are implemented in the `vec<Type, N>` (statically allocated) and `vec<Type>` (dynamically allocated) classes, while matrices are implemented in the `mat<Type, N, K>` and `mat<Type>` classes. `algebra_types.h` defines common matrix types to make them easily available: `vecX` and `matX` where X is 2,3,4 for real vectors and matrices, and the corresponding `cvecX` and `cmatX` for complex vectors and matrices.

```cpp
#include "algebra_types.h"
using namespace th;

int main() {
	
    // Initialize a 4-dimensional vector
    vec4 v = {1, 2, 3, 4};
    
    // Initialize a 4x4 matrix to the identity
    mat4 A = mat4::identity();
    
    // Initialize a 4-dimensional vector to all ones
    vec<> w = vec<real>(4, 1.0);
}
```

For dynamically allocated vectors and matrices, you may use the `vec<>` and `mat<>` for real elements and `cvec` and `cmat` for complex elements. Lower level linear algebra routines are implemented in the `algebra` namespace, for generic data structures.

### Vectors
Theoretica uses a hybrid paradigm between OOP and functional programming which tries to emulate **mathematical notation** when convenient. This means that generally mathematical structures have both **methods and pure functions** to operate on them.

```cpp
#include "vec.h"
using namespace th;


int main() {

    vec<> v1 = {1, 2, 3, 4};
    vec<> v2 = {3, 1, 2, 5};
    
    real a = v1 * v2;
    real b = dot(v1, v2);
    real c = v1.dot(v2);
    
    // ...
```

Here `dot()` and `operator*` do the same thing, they compute the dot product of the two vectors.

```cpp
    vec3 v = {1, 2, 3};
    vec3 w = {3, 1, 2};
    
    vec3 r1 = cross(v, w);
    vec3 r2 = v.cross(w);
```

The `cross()` function computes the cross product between the two 3-dimensional vectors, be it as a function or a method.

```cpp
    vec3 v = {1, 2, 3};
    vec3 q = {0, 5, 2};
    vec3 w = 3 * v + q;
```

Theoretica uses operator overloading to write equations for the vectors.

```cpp
    // Compute the norm and square norm of the vector
    real a = v.norm();
    real b = v.sqr_norm();
    
    // Compute the normalized vector
    vec3 z = v.normalized();
    
    // Normalize the vector itself
    v.normalize();
    algebra::normalize(v);
```

The `vectorized` namespace contains useful functions to compute mathematical functions on all of the elements of a vector, speeding up execution through parallelization.

```cpp
    // Compute the square root of all the elements in v
    w = vectorized::sqrt(v);
```

### Matrices

Matrices are implemented in the `mat<Type, N, K>`, where the case `N = 0, K = 0` is used for dynamically allocated matrices which may change size at runtime. The default constructor creates a matrix with all zero entries.

```cpp
    // Construct a 4x4 identity matrix
    mat4 A = mat4::identity();
    mat<> B = mat<real>::identity(4, 4);
    
    // Matrix-vector multiplication
    vec4 w = A * v;
```

The `mat` class implements common matrix operations and provides matrix initialization to different types of matrices of common use, such as rotation matrices in 2D and 3D:

```cpp
    // Rotation by PI/3 around an arbitrary axis
    mat3 R = mat3::rotation_3d(PI / 3, {1, 2, 3});
    
    // Diagonal matrix
    mat4 D = mat4::diagonal(2.0);
```

Like vectors, matrices have a wide range of methods and functions to operate on them:

```cpp
    // Rotation matrix around the x axis
    mat3 A = mat3::rotation_3d_xaxis(PI);
    
    // Matrix determinant
    real d1 = algebra::det(A);
    real d2 = A.det();
    
    // Matrix trace
    real t1 = algebra::trace(A);
    real t2 = A.trace();
```

## Complex numbers
Different complex algebras are implemented in the `complex` module. Complex numbers in algebraic form $a + ib$ are implemented in the `complex<Type>` class, while complex numbers in exponential form $\rho e^{i\theta}$ are implemented in the `phasor` class. Quaternions are available using the `quat<Type>` class.

```cpp
    // Initialize a complex number as 1 + i
    complex<> z = {1, 1};
    
    // Construct a complex rotor of angle PI / 2
    complex<> w = complex<>::rotor(PI / 2);
    
    // Complex multiplication and conjugation
    complex<> s = z * conjugate(w);
    
    // Convert a complex number to exponential form
    phasor phi = z;
    
    // Construct a quaternion 1 + 2i + 3j + 4k
    quat<> q = {1, 2, 3, 4};
```

The `complex_analysis.h` header contains common complex functions, such as the square root, trigonometric functions and so on (where the main branch is taken):

```cpp
    complex<> r = sqrt(z) + sin(w);
    real x = Re(r + inverse(z)) / Im(r);
```

Quaternions may be also used to rotate vectors in 3 dimensions:

```cpp
    // Return a quaternion representing a rotation
    // around an arbitrary axis
    quat<> r = quat<>::rotation(PI, v);

    // Rotate a vector around an arbitrary axis using quaternions
    vec3 w = quat<>::rotate(w, PI / 2.0, v);
```

## Derivatives and integrals
Theoretica provides several functions to approximate the derivatives (using finite differences) and integrals (using quadrature methods). These functions are defined inside the `calculus` module. Numerical methods for derivatives use the `deriv_` prefix, while numerical integrals use the `integral_` prefix.

```cpp
real f(real x) {
    return x * sqrt(x);
}

int main() {

    // Compute the derivative of f at x = 2.0
    // using a generic method.
    real Dx = deriv(f, 2.0);
    
    // Compute the second derivative of f at x = 2.0
    // using a generic method.
    real Dx2 = deriv2(f, 2.0);
    
    // Compute the integral of f in [3.0, 5.0]
    // using a generic method.
    real I = integral(f, 3.0, 5.0);
    
    // Compute the integral of f in [3.0, 5.0]
    // using Romberg's method.
    real I = integral_romberg(f, 3.0, 5.0);
}
```

Several specific methods are available for both derivatives and integrals. The forward, backward, central and Ridder's methods are used for derivatives, while trapezoid, Simpson's, Romberg's and different Gauss' methods are available for integrals. Derivatives may also be computed using **automatic differentiation**, which will be introduced in the next paragraph.

## Automatic Differentiation
Theoretica also implements automatic differentiation of univariate and multivariate real functions through the `autodiff` module. Univariate differentiation uses the `dual` class, while multivariate differentiation uses the `multidual<N>` class with typedefs `d_real<N>` and `d_vec<N>` for ease of use.

This examples illustrates using univariate differentiation:
```cpp
dual f(dual x) {
    return x * sqrt(x);
}

int main() {
	
    real Dx = f(dual(2.0, 1.0)).Dual();
}
```

Here the real of the dual number contains the value of $f(2.0)$ while the dual part contains the value of $f'(2.0)$.

The following example uses multivariate differentation:

```cpp

d_real<2> g(d_vec<2> v) {
    return v[0] * sqrt(v[1]);
}

int main() {
	
    vec2 v = {1, 2};
    
    // Compute the gradient of f at v
    vec2 w = gradient(f, v);
}
```

The gradient is automatically computed at the given point and returned.
Additional differential operators such as `jacobian()`, `divergence` and `curl` are also available in the `autodiff.h` header.


## Statistics
Theoretica provides routines for descriptive and inferential statistics inside the `statistics` module. Common estimators such as the sample mean and variance are defined for generic containers of data.

```cpp
    vec<> v = {1, 2, 4, 3, 1, 2};
    
    real m = mean(v);
    real var = sample_variance(v);
```

### Probability distributions
Many common probability distribution functions are defined in the `distributions.h` header, inside the `distribution` namespace.

```cpp
  real x1 = distribution::gaussian(1, 2.1, 0.7);
  real x2 = distribution::binomial(1, 3, 0.75);
  real x3 = distribution::poisson(2, 3);
```
### Linear regression
The `regression` namespace also contains routines for linear regression, in particular the `linear_model` class provides simple access to linear regression with or without uncertainties over the coordinates:

```cpp
    vec<> X = {1.1, 1.9, 2.7, 4.2};
    vec<> Y = {2.3, 3.8, 6.5, 7.6};
    
    // Compute the linear regression between X and Y
    // with uncertainty 0.2 over X
    auto l = regression::linear_model(X, Y, 0.2);
    
    // Get the computed p-value
    real P = l.p_value;
```

## Pseudorandom generators
Theoretica implements different pseudorandom generators made available through the `PRNG` class:

```cpp
    uint64_t seed = 123456789;
	
    // Xoshiro256++ generator
    PRNG g1 = PRNG::xoshiro(seed);
    
    // Wyrand generator
    PRNG g2 = PRNG::wyrand(seed);
    
    // Discard first 10'000 numbers
    g1.discard(10000);
    
    // Generate a random number
    uint64_t n = g1();
```

These PRNGs can be used through the `pdf_sampler` class to sample probability distribution functions:

```cpp
    // Standard Gaussian sampler
    pdf_sampler gauss = pdf_sampler::gaussian(g, 0, 1);
    
    // Sample a single value
    real x = gauss();
    
    // Fill a vector with random values
    vec<> v = vec<>(1000);
    gauss.fill(v);
```

## Polynomials
Polynomials may be easily manipulated using the `polynomial<Type>` class. The coefficients of the polynomial are passed in ascending degree:

```cpp
    // x + 2x^3
    polynomial<real> P = {0, 1, 0, 2};
    
    // Evaluate the polynomial at x = 1.0
    real r = P(1.0);
    
    // Divide the polynomial by x
    P /= polynomial({0, 1});
    
    polynomial<> P1 = {1, 2};
    polynomial<> P2 = {3, 4};

    // Linear combination of polynomials
    polynomial<> P3 = P1 + 2 * P2;
```

The order of the polynomial may be computed by using `.find_order()`, while a polynomial may be trimmed of zero higher order coefficients using `.trim()`.
