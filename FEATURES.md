# List of features
This is a non-comprehensive list of features. Elements marked with an 'x' have been fully implemented.

A functionality is considered "_fully implemented_" when it has at least one implementation which is architecture independent that passed all test cases and has no performance issues. Every function may have an **architecture independent** implementation, possibly of **arbitrary precision**, and one or more **architecture dependent** implementations, eventually making use of Assembly instructions.

- Real analysis (_real_analysis.h_)
  - [x] sqrt
  - [x] cbrt
  - [x] square
  - [x] cube
  - [x] abs
  - [x] sgn
  - [ ] max - Test cases
  - [ ] min - Test cases
  - [ ] clamp - Test cases
  - [ ] floor - Test cases
  - [ ] fract - Test cases
  - [ ] **log2** - **(only hardware implementation)**
  - [ ] log10 - **(only hardware implementation)**
  - [ ] ln - **(only hardware implementation)**
  - [x] pow - (may be improved)
  - [x] fact
  - [ ] **powf** - **(only hardware implementation)**
  - [ ] **exp** - **(only hardware implementation)**
  - [ ] sin - **(only hardware implementation)**
  - [ ] cos - **(only hardware implementation)**
  - [ ] tan - **(only hardware implementation)**
  - [ ] cot - **(only hardware implementation)**
  - [ ] **atan** - Low precision
  - [ ] asin - Low precision
  - [ ] acos - Low precision
  - [ ] **atan2** - Low precision
  - [ ] sinh - Depends on exp
  - [ ] cosh - Depends on exp
  - [ ] tanh - Depends on exp
  - [ ] coth - Depends on exp
  - [x] binomial_coeff
  - [x] radians
  - [x] degrees

- Vector algebra (_algebra/vector.h_) - **NO TEST CASES**
  - [x] Addition and subtraction (_operator+_, _-_)
  - [x] +=, -=, *= operators
  - [x] Scalar multiplication (_operator*_)
  - [x] Dot product (_operator*_)
  - [x] Cross product (_cross_)
  - [x] Norm (_lenght_)
  - [x] Normalization (_normalize_, _normalized_)

- Matrix algebra (_algebra/mat.h_) - **NO TEST CASES**
  - [x] Addition and subtraction (_operator+_, _-_)
  - [x] Scalar multiplication (_operator*_)
  - [x] Matrix-Vector product (_operator*_, _transform_)
  - [x] Matrix-Matrix product (_operator*_, _transform_)
  - [x] +=, -=, *= operators
  - [x] Transposition (_transpose_, _transposed_)
  - [x] Dot product of two vectors (_dot_)
  - [x] Matrix types (_is_square_, _is_diagonal_, _is_symmetric_)
  - [x] Transformation matrices (_identity_, _translation_, _rotation_4x4_, _rotation_3x3_, _scaling_)
  - [x] Perspective and Ortho matrices - Need testing
  - [x] Inversion
  - [x] Determinant of a generic matrix
  - [x] Determinant of 2x2 and 3x3 matrices

- Complex numbers (_complex/complex.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [x] Modulus
  - [x] Rotor construction
  - [ ] Argument - (depends on atan2 implementation)
  - [x] Conversion to other types (matrix form, to vector)

- Complex functions (_complex/complex_functions.h_) - **NO TEST CASES**
  - [x] square
  - [x] cube
  - [ ] exp - Depends on real function
  - [x] abs
  - [ ] sin - Depends on real function
  - [ ] cos - Depends on real function
  - [ ] tan - Depends on real function
  - [x] sqrt
  - [ ] ln - Depends on real function
  - [ ] asin - Depends on real function
  - [ ] acos - Depends on real function
  - [ ] atan - Depends on real function

- Phasors (_complex/phasor.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conversion from and to algebraic complex

- Quaternions (_complex/quat.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Norm
  - [x] Conjugate and inverse
  - [x] Vector transformation
  - [x] Matrix form
  - [x] Rotation construction

- Polynomials (_polynomial.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication
  - [x] Differentiation
  - [x] Integration
  - [ ] Division

- Calculus (_calculus/derivation.h_, _calculus/integration.h_) - **NO TEST CASES**
  - [x] Integral approximation using midpoint, trapezoid and Simpson
  - [ ] Romberg integral approximation
  - [ ] Monte Carlo integral approximation
  - [x] Derivative approximation

- Approximation of roots and extrema (_approx/roots.h_, _approx/extrema.h_) - **NO TEST CASES**
  - [x] Approximation of roots using Newton, bisection, Steffensen and Chebyshev
  - [x] Approximation of extrema using Newton, bisection and golden section search

- Dual numbers (_autodiff/dual.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [x] Conversion from and to vector
  - [x] Matrix form

- Multidual numbers (_autodiff/multidual.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [ ] Conversion from and to vector
  - [ ] Matrix form

- Dual functions and Automatic Differentiation (_autodiff/dual_functions.h_) - **NO TEST CASES**
  - [x] square
  - [x] cube
  - [x] pow
  - [x] sqrt
  - [ ] sin - Depends on real function
  - [ ] cos - Depends on real function
  - [ ] tan - Depends on real function
  - [ ] cot - Depends on real function
  - [ ] exp - Depends on real function
  - [ ] ln - Depends on real function
  - [ ] log2 - Depends on real function
  - [ ] log10 - Depends on real function
  - [ ] asin - Depends on real function
  - [ ] acos - Depends on real function
  - [ ] atan - Depends on real function
  - [x] abs

- Interpolation (_interpolation/polyn_interp.h_, _interpolation/spline_interp.h_) - **NO TEST CASES**
  - [x] Polynomial interpolation
  - [x] Linear interpolation
  - [x] Inverse linear interpolation
  - [x] Remapping
  - [x] Quadratic Bezier
  - [x] Cubic Bezier
  - [x] Generic Bezier

- Pseudorandom numbers (_pseudorandom/pseudorandom_algo.h_, _pseudorandom/prng.h_) - **NO TEST CASES**
  - [x] PRNG class
  - [x] Random real numbers on a range
  - [x] Congruential generator
  - [ ] Mersenne twister
  - [ ] Xorshift
  - [ ] Lagged Fibonacci

- Statistical functions (_statistics/statistics.h_) - **NO TEST CASES**
  - [x] Mean and weighted mean
  - [x] Sum and product error propagation
  - [x] Total sum of squares
  - [x] Variance
  - [x] Covariance
  - [x] Standard deviation
  - [x] Standard deviation of the mean
  - [x] Correlation coefficient
  - [x] Least squares linearization
  - [x] Chi square of a linearization

- Distributions (_statistics/distributions.h_) - **NO TEST CASES**
  - [x] Likelihood and log likelihood
  - [x] Gaussian
  - [x] Bernoulli
  - [x] Poisson
  - [x] Binomial
  - [x] Log-normal
  - [x] Exponential
  - [x] Cauchy
  - [x] Breit Wigner

- Random number-based algorithms (_pseudorandom/randstat.h_) - **NO TEST CASES**
  - [x] Try-and-catch number generation following a distribution
  - [ ] Inverse cumulative function method
  - [ ] Metropolis-Hastings Monte Carlo Markov Chain method
