# List of features
This is a non comprenhensive list of features. Elements marked with an 'x' have been fully implemented.

A functionality is considered "_fully implemented_" when it has at least one implementation which is architecture independent that passed all test cases and has no performance issues. Every function may have an **architecture independent** implementation, possibly of **arbitrary precision**, and one or more **architecture dependent** implementations, eventually making use of Assembly instructions.

- Real analysis (_real_analysis.h_)
  - [ ] **sqrt**
  - [ ] cbrt
  - [x] square
  - [x] cube
  - [x] abs
  - [x] sgn
  - [x] max
  - [x] min
  - [x] clamp
  - [ ] fyl2x - **(MSVC)**
  - [ ] f2xm1 - **(MSVC)**
  - [ ] **log2** - **(only hardware implementation)**
  - [ ] log10 - **(only hardware implementation)**
  - [ ] ln - **(only hardware implementation)**
  - [x] pow - (may be improved)
  - [x] fact
  - [ ] exp_approx - (needs a better implementation)
  - [ ] **powf_approx** - (needs a better implementation)
  - [ ] **exp**
  - [ ] sin - **(only hardware implementation)**
  - [ ] cos - **(only hardware implementation)**
  - [ ] tan - **(only hardware implementation)**
  - [ ] cot - **(only hardware implementation)**
  - [ ] **atan**
  - [ ] asin
  - [ ] acos
  - [ ] **atan2**
  - [ ] sinh
  - [ ] cosh
  - [ ] tanh
  - [ ] coth
  - [ ] binomial_coeff - (lacks test cases)
  - [x] radians
  - [x] degrees

- Vector algebra (_algebra/vector.h_) - **NO TEST CASES**
  - [x] Addition and subtraction (_operator+_, _-_)
  - [ ] +=, -=, *= operators
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
  - [ ] +=, -=, *= operators
  - [x] Transposition (_transpose_, _transposed_)
  - [x] Dot product of two vectors (_dot_)
  - [x] Matrix types (_is_square_, _is_diagonal_, _is_symmetric_)
  - [x] Transformation matrices (_identity_, _translation_, _rotation_4x4_, _rotation_3x3_, _scaling_)
  - [ ] Perspective and Ortho matrices - Need testing
