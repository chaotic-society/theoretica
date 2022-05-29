
///
/// @file theoretica.h General include file
///
/// This file will include core header files.
/// If `THEORETICA_INCLUDE_ALL` is defined, it will include all header files of the library
///

#ifndef THEORETICA_H
#define THEORETICA_H

// Core headers
#include "./core/constants.h"
#include "./core/error.h"

// Real functions
#include "./core/real_analysis.h"

// Vector and matrix algebra
#include "./algebra/vec.h"
#include "./algebra/mat.h"

// Complex and quaternion classes
#include "./complex/complex.h"
#include "./complex/quat.h"


// All other headers
#ifdef THEORETICA_INCLUDE_ALL

// Ratio
#include "./core/ratio.h"

// Bit operations
#include "./core/bit_op.h"

// Data sets
#include "./core/vec_buff.h"

// Interpolation
#include "./interpolation/spline_interp.h"
#include "./interpolation/polyn_interp.h"
	
// Statistical functions
#include "./statistics/statistics.h"
#include "./statistics/distributions.h"

// Roots and extrema approximation of real functions
#include "./approx/roots.h"
#include "./approx/extrema.h"

// Complex and quaternion functions
#include "./complex/complex_analysis.h"
#include "./complex/phasor.h"

// Derivative and integral approximation
#include "./calculus/derivation.h"
#include "./calculus/integration.h"
#include "./calculus/taylor_expansion.h"

// Polynomial class
#include "./polynomial/polynomial.h"

// Dual numbers and automatic differentiation
#include "./autodiff/dual.h"
#include "./autodiff/dual_functions.h"
#include "./autodiff/multidual.h"
#include "./autodiff/multidual_functions.h"
#include "./autodiff/dual2.h"
#include "./autodiff/dual2_functions.h"
#include "./autodiff/autodiff.h"

// Pseudorandom number generation
#include "./pseudorandom/pseudorandom_algo.h"
#include "./pseudorandom/prng.h"
#include "./pseudorandom/rand_dist.h"

// Monte Carlo methods
#include "./pseudorandom/montecarlo.h"

#endif

/// @namespace th Alias for the theoretica namespace
#ifndef THEORETICA_NO_NAMESPACE
namespace th = theoretica;
#endif

#endif
