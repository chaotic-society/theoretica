<img alt="GitHub code size in bytes" src="https://img.shields.io/github/languages/code-size/mattiaisgro/uroboro"> <img alt="GitHub issues" src="https://img.shields.io/github/issues/mattiaisgro/uroboro">
# Uroboro
A general purpose mathematical header-only library in C++ intended for graphical and physical simulations. Includes approximations of common math functions with Assembly enhancements, vector and matrix calculations, as well as complex numbers, quaternions and statistical functions (see [Functionalities](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Functionalities) for more).

Matrices are implemented column-first, for direct use in OpenGL.
See [FEATURES.md](https://github.com/mattiaisgro/uroboro/blob/master/FEATURES.md) for a list of the features of the library and [PROJECT_STRUCTURE.md](https://github.com/mattiaisgro/uroboro/blob/master/PROJECT_STRUCTURE.md) for an overview of the project file structure.

## How to use
Include `uroboro.h` to be able to use common real functions, the `vec<N>` and `mat<N, K>` classes.
All functions and classes are defined inside the `uroboro` namespace (and the `umath` alias if `UROBORO_NO_NAMESPACE` is _not_ defined).

Define `UROBORO_INCLUDE_ALL` before including `uroboro.h` to automatically include all header files, otherwise only basic functionality will be included (vectors, matrices, complex numbers and common functions).

## Functionalities
- Common math operations (trig. functions, sqrt, exp, ...)
- Vector operations (N-dimensional)
- Matrix operations (with transformations, including rotations)
- Complex numbers (with common complex functions), Quaternions and Phasors
- Wide range of statistical functions (including least squares linearization)
- Approximation of roots and extrema
- Derivative and integral approximations (including Runge-Kutta of 4th order)
- Interpolation between points (Linear interpolation, Bezier)

## Examples
The file `example.cpp` gives example usage of the library.

### Vector usage
Declare a vector using the `vec<N>` type or predefined types `vec2`, `vec3` or `vec4`.

```c

using namespace umath; // same as uroboro

vec3 v = {1, 2, 3};
vec3 u = {1, 0, 1};

real d = v * u; // Dot product
vec w = v.cross(u); // Cross product

// 12 dimensional vector
vec<12> a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};

real b = a.lenght(); // Lenght of a vector (same as magnitude())
vec3 a_n = a.normalized(); // Return the normalized vector
a.normalize(); // Normalize the vector itself

```

### Matrix usage
Declare a matrix using the `mat<N, K>` type (where `N` is the number of columns and `K` the number of rows) or the predefined square matrices types 
`mat2`, `mat3` or `mat4`.

```c
mat4 A = mat4::identity(); // A is a 4x4 identity matrix
mat3 B = mat3::diagonal(3.0); // B is a 3x3 diagonal matrix
mat3 C = mat4::rotation_3x3(PI, u); // C is a 3x3 rotation matrix around the given axis

vec3 l = B * v;
vec3 t = C * w;
```

### Complex number usage
Declare a complex number using the `complex` class.
```c
complex z1 = {1, 1};
complex z2 = complex::rotor(TAU); // A complex rotor

complex z3 = z1 * z2; // Complex multiplication
complex z4 = z1.conjugate(); // Complex conjugate
complex z5 = ln(z3) + atan(z3);
```
Common complex functions such as `sqrt`, `ln`, complex trigonometric functions and more are also supported (`complex_analysis.h`).

### Quaternion usage
Declare a quaternion using the `quat` class.
```c
quat q = {1, 2, 3, 4};

quat p = q.inverse(); // Return the inverse of a quaternion
quat p_n = p.normalized(); // Return the normalized quaternion

// Return a quaternion representing a rotation around an arbitrary axis
quat r = quat::rotation(PI, v);

// Rotate a vector around an arbitrary axis using quaternions
vec3 v_1 = quat::rotate(w, PI / 2.0, v);
```

### Polynomials
Polynomial storage and manipolation is implemented in `polynomial.h`. Derivation and integration of polynomials is also supported (`derivation.h`, `integration.h`). The header file `approx.h` defines many functions to approximate the roots of arbitrary functions and polynomials.
```c
polynomial P1 = {1, 1, -1};
polynomial P2 = {2, 3, 6, 1};
polynomial P3 = P1 * P2 + differentiate_polynomial(P1);
print_polynomial(P3);
```

### Statistics usage
The header `statistics/statistics.h` defines many statistical functions of common use. The `vec_buff` type (alias for `std::vector<real>`) is used to store data.
```c
vec_buff data = {1.0, 1.1, 1.3, 0.9, 1.0, 0.9, 0.8};

real x_m = mean(data);
real sigma = smpl_stdev(data); // Standard deviation of a sample (same as sample_standard_deviation)
real sigma_m = smpl_stdom(data); // STDOM of a sample (same as sample_mean_standard_deviation)

vec_buff X = {0.9, 2.1, 3.2, 3.9, 5.2};
vec_buff Y = {0.11, 0.22, 0.29, 0.41, 0.53};
vec_buff W = {0.01, 0.01, 0.02, 0.01, 0.03}; // Weights of the errors on Y

real c = sample_covariance(X, Y);

// Linearization using least squares
real intercept = lst_sqrs_lin_intercept(X, Y);
real slope = lst_sqrs_lin_slope(X, Y);

real intercept2 = lst_sqrs_weight_lin_intercept(X, Y);
real slope2 = lst_sqrs_weight_lin_slope(X, Y);

real r = sample_correlation_coefficient(X, Y);

```

Many probability distribution functions are implemented in `statistics/distributions.h` and can be used through the `distribution` namespace. These include Gaussian, Binomial, Log-Normal, Poisson, Bernoulli, Cauchy and Breit Wigner.
```c
real x1 = distribution::gaussian(1, 2.1, 0.7);
real x2 = distribution::binomial(1, 3, 0.75);
real x3 = distribution::poisson(2, 3);
```

## Macros
These are the macros that can be defined to change the library's behaviour:
- **UROBORO_INCLUDE_ALL** - Including `uroboro.h` will include _all_ header files
- **UROBORO_BRANCHLESS** - Branchless implementations will be preferred when possible (and computationally convenient)
- **UROBORO_FORCE_BRANCHLESS** - Branchless implementations will always be used whenever possible
- **UROBORO_THROW_EXCEPTIONS** - Exceptions will be thrown and errno set on error (by default errno is set and NaN is returned)
- **UROBORO_ONLY_EXCEPTIONS** - Exceptions will be thrown on error (without modifying errno)
- **UROBORO_NO_NAMESPACE** - The `umath` namespace alias will not be defined
- **UROBORO_X86** - Assembly x86 implementations will be used whenever possible
- **UROBORO_FLOAT_PREC** - Floating point precision (`float`) will be used for the `real` type (by default `double` is used)
- **UROBORO_LONG_DOUBLE_PREC** - Long double precision (`long double`) will be used
- **UROBORO_ARBITRARY_PREC** - Arbitrary precision will be used
- See `constants.h` for more specific defines

