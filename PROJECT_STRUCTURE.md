# Project structure

### The _src_ folder
Where source code is stored.

- _src/algebra_ - Linear algebra classes
  - _mat.h_ - Matrix class
  - _vec.h_ - Vector class
  
- _src/calculus_ - Integral and differential calculus
  - _derivation.h_ - Derivative approximation
  - _integration.h_ - Integral approximation
  
- _src/complex_ - Complex numbers
  - _complex.h_ - Complex number class in algebraic form
  - _complex_analysis.h_ - Complex functions
  - _phasor.h_ - Complex number class in exponential form
  - _quat.h_ - Quaternion class

- _src/autodiff_ - Automatic Differentiation using dual numbers

- _src/other_ - Miscellaneous files
  - _interpolation.h_ - Interpolation functions on vectors

- _src/statistics_ - Statistical functions
  - _distributions.h_ - Probability distribution functions
  - _statistics.h_ - Common statistical functions

- _approx/roots.h_ - Approximation of roots of real functions
- _approx/extrema.h_ - Approximation of extrema of real functions
- _constants.h_ - Constants header
- _example.cpp_ - Example file showcasing common use functions
- _function.h_ - Function types header
- _polynomial.h_ - Polynomial class
- _real_analysis.h_ - Common real functions
- _uroboro.h_ - Header file including headers of common use
- _utility.h_ - Printing functions
- _vec_buff.h_ - Dynamic vector class

## The _test_ folder
Where test cases are stored.

- _test_*.cpp_ - Test cases files
- _test.h_ - Testing library code
