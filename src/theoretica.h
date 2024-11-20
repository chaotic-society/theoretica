
///
/// @file theoretica.h General include file
///
/// If `THEORETICA_INCLUDE_BASE` is defined, only core headers will be included
///

#ifndef THEORETICA_H
#define THEORETICA_H

// Core headers
#include "core/constants.h"
#include "core/error.h"
#include "core/core_traits.h"

// Real functions
#include "core/real_analysis.h"

// Vector and matrix algebra
#include "algebra/algebra.h"
#include "algebra/algebra_types.h"
#include "algebra/transform.h"
#include "algebra/vec.h"
#include "algebra/mat.h"
#include "algebra/distance.h"

// Complex and quaternion classes
#include "complex/complex.h"
#include "complex/quat.h"


// All other headers
#ifndef THEORETICA_INCLUDE_BASE

// Special functions
#include "core/special.h"

// Vectorized functions
#include "algebra/parallel.h"

// Ratio
#include "core/ratio.h"

// Bit operations
#include "core/bit_op.h"

// Data sets
#include "core/dataset.h"

// Interpolation
#include "interpolation/splines.h"
#include "interpolation/polynomial.h"
	
// Statistical functions
#include "statistics/statistics.h"
#include "statistics/distributions.h"
#include "statistics/errorprop.h"
#include "statistics/histogram.h"
#include "statistics/regression.h"

// Roots and extrema approximation of real functions
#include "optimization/roots.h"
#include "optimization/extrema.h"
#include "optimization/multi_extrema.h"
#include "optimization/multi_roots.h"

// Complex and quaternion functions
#include "complex/complex_analysis.h"
#include "complex/phasor.h"
#include "complex/complex_types.h"

// Derivative and integral approximation
#include "calculus/deriv.h"
#include "calculus/integral.h"
#include "calculus/ode.h"
#include "calculus/taylor.h"

// Polynomial class
#include "polynomial/polynomial.h"
#include "polynomial/orthogonal.h"

// Dual numbers and automatic differentiation
#include "autodiff/dual.h"
#include "autodiff/dual_functions.h"
#include "autodiff/multidual.h"
#include "autodiff/multidual_functions.h"
#include "autodiff/dual2.h"
#include "autodiff/dual2_functions.h"
#include "autodiff/autodiff.h"

// Pseudorandom number generation
#include "pseudorandom/pseudorandom.h"
#include "pseudorandom/prng.h"
#include "pseudorandom/sampling.h"

// Monte Carlo methods
#include "pseudorandom/montecarlo.h"

// Fast Fourier transform
#include "signal/fft.h"

#endif

#endif
