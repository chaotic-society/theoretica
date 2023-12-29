# List of features
This is a non-comprehensive list of features. Elements marked with an 'x' have been fully implemented.

A functionality is considered "fully implemented" when it has at least one implementation which is architecture independent that passed all test cases and has no performance issues. Every function may have an **architecture independent** implementation, possibly of **arbitrary precision**, and one or more **architecture dependent** implementations, eventually making use of Assembly instructions.

- Real analysis (_core/real_analysis.h_)
  - [x] sqrt
  - [x] cbrt
  - [x] square
  - [x] cube
  - [x] isqrt
  - [x] icbrt
  - [x] abs
  - [x] sgn
  - [x] max
  - [x] min
  - [x] clamp
  - [x] floor
  - [x] fract
  - [x] log2
  - [x] log10
  - [x] ln
  - [x] ilog2
  - [x] pow
  - [x] powf
  - [x] fact
  - [x] double_fact
  - [x] exp
  - [x] sin
  - [x] cos
  - [x] tan
  - [x] cot
  - [x] atan
  - [x] asin
  - [x] acos
  - [x] atan2
  - [x] sinh
  - [x] cosh
  - [x] tanh
  - [x] coth
  - [x] binomial_coeff
  - [x] radians
  - [x] degrees
 
- Rational class (_core/ratio.h_)
  - [x] Addition
  - [x] Multiplication
  - [x] Division
  - [x] Evaluation
  - [x] Conversion to string

- Vector algebra (_algebra/vector.h_)
  - [x] Addition and subtraction (_operator+_, _-_)
  - [x] +=, -=, *= operators
  - [x] Scalar multiplication (_operator*_)
  - [x] Dot product (_operator*_)
  - [x] Cross product (_cross_)
  - [x] Norm (_lenght_)
  - [x] Normalization (_normalize_, _normalized_)
  - [x] Euclidean distance
  - [x] Lp-norm

- Matrix algebra (_algebra/mat.h_, _algebra/algebra.h_)
  - [x] Addition and subtraction (_operator+_, _-_)
  - [x] Scalar multiplication (_operator*_)
  - [x] Matrix-Vector product (_operator*_, _transform_)
  - [x] Matrix-Matrix product (_operator*_, _transform_)
  - [x] +=, -=, *= operators
  - [x] Transposition (_transpose_, _transposed_)
  - [x] Dot product of two vectors (_dot_)
  - [x] Matrix types (_is_square_, _is_diagonal_, _is_symmetric_)
  - [x] Transformation matrices (_identity_, _translation_, _rotation_4x4_, _rotation_3x3_, _scaling_)
  - [x] Perspective and Ortho matrices
  - [x] Inverse of a generic matrix
  - [x] Determinant of a generic matrix

- Distances and norms (_algebra/distance.h_)
  - [x] Lp norm
  - [x] L1, L2, Linf norm
  - [x] Euclidean distance
  - [x] Minkowski distance
  - [x] Hermitian distance
  - [x] Manhattan distance
  - [x] Chebyshev distance
  - [x] Discrete distance
  - [x] Hamming distance

- Complex numbers (_complex/complex.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [x] Modulus
  - [x] Rotor construction
  - [x] Argument
  - [x] Conversion to other types (matrix form, to vector)

- Complex functions (_complex/complex_functions.h_) - **NO TEST CASES**
  - [x] square
  - [x] cube
  - [x] exp
  - [x] abs
  - [x] sin
  - [x] cos
  - [x] tan
  - [x] sqrt
  - [x] ln
  - [x] asin
  - [x] acos
  - [x] atan

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

- Polynomials (_polynomial/polynomial.h_)
  - [x] Addition and subtraction
  - [x] Multiplication
  - [x] Differentiation
  - [x] Integration
  - [x] Division

- Orthogonal polynomials (_polynomial/ortho_polyn.h_) - **NO TEST CASES**
  - [x] Laguerre polynomials
  - [x] Legendre polynomials
  - [x] Hermite polynomials
  - [x] Chebyshev polynomials

- Calculus (_calculus/derivation.h_, _calculus/integration.h_, _calculus/taylor.h_)
  - [x] Integral approximation using midpoint, trapezoid and Simpson
  - [x] Romberg integral approximation
  - [x] Monte Carlo integral approximation (Crude, Hit-or-Miss)
  - [x] Derivative approximation using central, forward, backward methods
  - [x] Ridder's derivative approximation of arbitrary order
  - [x] Second derivative approximation
  - [x] Ordinary Differential Equations numerical integration methods
  - [x] Linear and quadratic Taylor series expansions of generic functions

- Approximation of roots and extrema (_optimization/roots.h_, _optimization/extrema.h_) - **NO TEST CASES**
  - [x] Approximation of roots using Newton, bisection, Steffensen and Chebyshev
  - [x] Approximation of extrema using Newton, bisection and golden section search

- Approximation of roots and extrema of multivariate functions (_optimization/multi_extrema.h_)
  - [x] Crude gradient descent
  - [x] Gradient descent with linear search

- Dual numbers (_autodiff/dual.h_, _autodiff/dual2.h_)
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [x] Conversion from and to vector
  - [x] Matrix form (only first order duals)

- Multidual numbers (_autodiff/multidual.h_) - **NO TEST CASES**
  - [x] Addition and subtraction
  - [x] Multiplication and division
  - [x] Conjugate and inverse
  - [ ] Conversion from and to vector
  - [ ] Matrix form

- Dual functions and Automatic Differentiation (_autodiff/dual_functions.h_, _autodiff/dual2_functions.h_)
  - [x] square
  - [x] cube
  - [x] pow
  - [x] sqrt
  - [x] sin
  - [x] cos
  - [x] tan
  - [x] cot
  - [x] exp
  - [x] ln
  - [x] log2
  - [x] log10
  - [x] asin
  - [x] acos
  - [x] atan
  - [x] abs

- Differential operators using automatic differentiation (_autodiff/autodiff.h_) - **NO TEST CASES**
  - [x] Gradient
  - [x] Directional derivative
  - [x] Jacobian
  - [x] Divergence
  - [x] Curl
  - [x] Laplacian
  - [ ] Hessian
  - [x] Sturm-Liouville

- Interpolation (_interpolation/polyn_interp.h_, _interpolation/spline_interp.h_) - **NO TEST CASES**
  - [x] Polynomial interpolation
  - [x] Linear interpolation
  - [x] Inverse linear interpolation
  - [x] Spherical interpolation
  - [x] Remapping
  - [x] Quadratic Bezier
  - [x] Cubic Bezier
  - [x] Generic Bezier
  - [x] Natural Cubic Splines

- Pseudorandom numbers (_pseudorandom/pseudorandom.h_, _pseudorandom/prng.h_)
  - [x] PRNG class
  - [x] Congruential generator
  - [x] Xoroshiro256++
  - [x] Splitmix64
  - [x] Wyrand

- Quasirandom numbers (_pseudorandom/quasirandom.h_) - **NO TEST CASES**
  - [x] Weyl sequence
  - [ ] Sobol sequence
  - [ ] Halton sequence

- Statistical functions (_statistics/statistics.h_)
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
  - [x] p-value

- Linear model regression (_statistics/regression.h_)
  - [x] Linear Ordinary Least Squares
  - [x] Linear Weighted Least Squares
  - [x] Linear Ordinary Least Squares through origin
 
- Error propagation (_statistics/errorprop.h_) - **NO TEST CASES**
  - [x] Automatic error propagation
  - [x] Covariance Matrix estimation
  - [x] Monte Carlo error propagation

- Distributions (_statistics/distributions.h_) - **NO TEST CASES**
  - [x] Likelihood and log likelihood
  - [x] Gaussian
  - [x] Bernoulli
  - [x] Poisson
  - [x] Binomial
  - [x] Multinomial
  - [x] Log-normal
  - [x] Exponential
  - [x] Cauchy
  - [x] Breit Wigner
  - [x] Laplace
  - [x] Erlang
  - [x] Pareto
  - [x] Chi-squared
  - [x] Student's t

- Histogram class (_statistics/histogram.h_) - **NO TEST CASES**
  - [x] Histogram initialization from data
  - [x] Frequency histogram to string format
  - [x] Running statistics

- Random number-based algorithms (_pseudorandom/rand_dist.h_, _pseudorandom/montecarlo.h_) - **NO TEST CASES**
  - [x] Try-and-catch number generation following a distribution
  - [x] pdf_sampler class
  - [x] Uniform sampling
  - [x] Gaussian sampling (Polar, Box-Muller, CLT)
  - [x] Exponential sampling
  - [x] Cauchy sampling
  - [x] Laplace sampling
  - [x] Pareto sampling
  - [x] Crude Monte Carlo
  - [x] Hit-or-miss Monte Carlo
  - [x] Metropolis-Hastings Monte Carlo Markov Chain method

- Special functions (_core/special.h_)
  - [x] Gamma function
  - [x] Pi function
  - [x] Beta function

- Operations on bits (_core/bit_op.h_) - **NO TEST CASES*
  - [x] 128-bit integer multiplication
  - [x] MUM bit mixing
  - [x] Bit rotation
