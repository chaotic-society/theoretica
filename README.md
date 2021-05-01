# Uroboro
A header-only math library for videogames, simulations and graphics in C++, with x86 Assembly enhancements.
Includes approximations of common math functions, vector and matrix calculations (see [Functionalities](https://github.com/mattiaisgro/uroboro/blob/master/README.md#Functionalities) for more).

Matrices are column-first, for direct use in OpenGL.

## How to use
Include `uroboro.h` and you will be able to use `vec<N>` and `mat<N, K>` classes.
All functions and classes are defined inside the `uroboro` namespace (and the `umath` alias if `UROBORO_NO_NAMESPACE` is _not_ defined).

Define `UROBORO_INCLUDE_ALL` before including `uroboro.h` and all header files will be included automatically, otherwise only basic functionality will be included.

## Functionalities
- Common math operations (trig. functions, sqrt, exp, ...)
- Vector operations
- Matrix operations (with transformations including rotations)
- Complex numbers and Quaternions
- Wide range of statistical functions (including least squares linearization)
- Derivative and integration approximations (including Runge-Kutta of 4th order)

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
```

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

### Statistics usage
The header `statistics` defines many statistical functions of common use. The `vec_buff` type (alias for `std::vector<real>`) is used to store data.
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
