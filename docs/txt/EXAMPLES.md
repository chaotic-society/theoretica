# Examples

## Including library headers

```c
#include "theoretica.h"

using namespace theoretica;
```

## Using vectors
Vector operations can be used through the `vec<N, T = real>` class, where `N` specifies the number of rows and `T` specifies the element type (`real` is used as default). Common vector types are defined as `vec2`, `vec3` and `vec4`.

Declare a vector:
```c
// vec3 is the same as vec<3, real>
vec3 v1 = {1, 2, 3};

// You can initialize a vector from a bracketed list
vec3 v2 = {0, 1, 0};
```

Basic vector operations:
```c
// Compute the length of a vector
real l = v1.length();

// Access the i-th element (by reference)

// You can use the vec::at() function
v1.at(1) = 2;

// Or the [i] operator
v1[1] = 2;

// Get and set the i-th element
l = v1.get(1);
v1.set(1, 2);

// Multiply a vector by a scalar
v1 = v1 * 2.0;
v2 = v2 / 3.0;

// Normalize a vector

// normalized() returns the normalized vector
// without changing the initial vector
vec3 v1_n = v1.normalized();

// Meanwhile normalize() normalizes the vector itself
v1.normalize()
```

Dot product:
```c
// You can use the class method
real d = v1.dot(v2);

// Or the two-argument function
d = dot(v1, v2);

// Or simply the * operator
d = v1 * v2;
```

Cross product:
```c
// Same thing with the cross product
vec3 c;
c = v1.cross(v2);
c = cross(v1, v2);

// But there is no cross product operator!
```

## Using matrices
Matrices are implemented in the `mat<N, K>` class, where N is the number of rows and K is the number of columns (colexicographical and column-first implementation). Common square matrix types are defines as `mat2`, `mat3`, `mat4`

Initialize a matrix:
```c
// 4x4 identity matrix
mat4 A1 = mat4::identity();

// Diagonal matrix with 3 on the diagonal
mat4 A2 = mat4(3.0);

// Zero matrix
mat3 B1 = mat3();
```

Basic matrix operations:
```c
// Access an element
B1.at(1, 1) = 2.0;

// Get or set and element
real a12 = A2.get(1, 2);
A1.set(a12, 1, 2);

// Scalar multiplication
A1 = (A2 * 2.0) / 3.0;
```

Vector and matrix multiplication:
```c
// Matrix-Vector multiplication

// You can use transform()
vec3 u1;
u1 = B1.transform(v1);

// Or simply the * operator
u1 = B1 * v1;

// Same with Matrix-Matrix multiplication
mat4 R1;
R1 = A1.transform(A2);
R1 = A1 * A2;
```

Matrix transposition:
```c
// transposed() returns the transposed matrix
// without modifying the original matrix
mat3 R2 = R1.transposed();

// Meanwhile transpose() transposes the matrix itself
R1.transpose();
```

Matrix inversion:
```c
// Compute the inverse of a matrix
mat4 I1 = A1.inverse();
```

Matrix determinant:
```c
// Compute the determinant
real D = I1.det();
```

Transformations:
```c
// Get the 3x3 rotation matrix around an arbitrary axis
mat3 Rot1 = mat3::rotation_3x3(PI2, vec3({1, 2, 3}));

// Rotation around the x axis (same for y and z)
mat3 Rot2 = mat3::rotation_x_3x3(PI4);

// 4x4 rotation matrix around the z axis of 10 degrees
// (radians(real) converts degrees to radians, degrees does the opposite)
mat4 Rot3 = mat4::rotation_z_4x4(radians(10));

// 4x4 translation matrix
mat4 T1 = mat4::translation(1, 2, 3);

// perspective(), ortho(), scaling_3x3(), scaling_4x4(),
// lookAt() functions are also supported
```

## Using complex numbers
Declare a complex number using the `complex` class.
```c
// You can initialize a complex number
// using a bracketed list
complex z1 = {1, 1};

// And you can also access its elements directly
z1.a = 1;
z2.b = 1;

// Or get its real and imaginary part
real x = z1.Re() + z1.Im();

// A complex rotor
complex z2 = complex::rotor(TAU);

// Complex multiplication
complex z3 = z1 * z2;

// Complex conjugate
complex z4 = z1.conjugate();

// Complex analysis functions
complex z5 = ln(z3) + atan(z3);
```
Common complex functions such as `sqrt`, `ln`, complex trigonometric functions and many more are also supported (`complex/complex_analysis.h`).

## Using quaternions
Declare a quaternion using the `quat` class.
```c
// Initialize a quaternion as 1 + 2i + 3j + 4k
quat q = {1, 2, 3, 4};

// Compute the inverse of a quaternion
quat p = q.inverse();

// Compute the normalized quaternion
quat p_n = p.normalized();

// Normalize the quaternion itself
p.normalize();

// Return a quaternion representing a rotation around an arbitrary axis
quat r = quat::rotation(PI, v);

// Rotate a vector around an arbitrary axis using quaternions
vec3 v_1 = quat::rotate(w, PI / 2.0, v);
```

## Using polynomials
Polynomials are stored and manipulated using the `polynomial<T = real>` class, where `T` is the type of the coefficients of the polynomial (`real` by default).
```c
// Represents the polynomial 1 + 2x + 3x^2
polynomial<real> P1 = {1, 2, 3};

// x + 2x^3
polynomial<real> P2 = {0, 1, 0, 2};

// Access the polynomial's coefficients
P2[1] = 2;

// Polynomial operations
polynomial<real> Pr = (P1 * P2) + P1;

// Evaluate the polynomial for a given value

// You can use the eval function
real res;
res = Pr.eval(1.0);

// Or simply the () operator
res = Pr(1.0);

// Find the (effective) order of a polynomial
real order = Pr.find_order();

// Trim unneeded higher-order null coefficients
Pr.trim();
```

## Using statistical functions
The header `statistics/statistics.h` defines many statistical functions of common use. The `vec_buff` type (alias for `std::vector<real>`) is used to store data. Many aliases are defined to shorten function names (e.g. `sample_standard_deviation` = `smpl_stdev`).

```c
vec_buff data = {1.0, 1.1, 1.3, 0.9, 1.0, 0.9, 0.8};

real x_m = mean(data);

// Standard deviation of a sample (same as sample_standard_deviation)
real sigma = smpl_stdev(data);

// STDOM of a sample (same as sample_mean_standard_deviation)
real sigma_m = smpl_stdom(data);

vec_buff X = {0.9, 2.1, 3.2, 3.9, 5.2};
vec_buff Y = {0.11, 0.22, 0.29, 0.41, 0.53};

// Weights of the errors on Y
vec_buff W = {0.01, 0.01, 0.02, 0.01, 0.03};

// Covariance of two datasets
real c = sample_covariance(X, Y);

// Linearization using least squares
real intercept = lst_sqrs_lin_intercept(X, Y);
real slope = lst_sqrs_lin_slope(X, Y);

real intercept2 = lst_sqrs_weight_lin_intercept(X, Y);
real slope2 = lst_sqrs_weight_lin_slope(X, Y);

// Correlation coefficient of two datasets
real r = sample_correlation_coefficient(X, Y);

```

Many probability distribution functions are implemented in `statistics/distributions.h` and can be used through the `distribution` namespace. These include Gaussian, Binomial, Poisson, Bernoulli, Cauchy, Chi-squared and Student's t, among others.
```c
real x1 = distribution::gaussian(1, 2.1, 0.7);
real x2 = distribution::binomial(1, 3, 0.75);
real x3 = distribution::poisson(2, 3);
```

The `pseudorandom/rand_dist.h` provides many algorithms to sample from different distributions, and the `pdf_sampler` class may be used to construct an object which automatically samples from a given distribution.
```c
// Construct a standard PRNG which uses the Wyrand algorithm
PRNG g = PRNG::wyrand();

// Sample a random value from a Gaussian distribution with mean 0 and sigma 1
real x = rand_gaussian(0, 1, g);
```
