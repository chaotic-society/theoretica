
///
/// @file theoretica.h General include file
///
/// If `THEORETICA_INCLUDE_BASE` is defined, only core headers will be included
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
#include "./algebra/distance.h"
#include "./algebra/geometry.h"

// Complex and quaternion classes
#include "./complex/complex.h"
#include "./complex/quat.h"


// All other headers
#ifndef THEORETICA_INCLUDE_BASE

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
#include "./optimization/roots.h"
#include "./optimization/extrema.h"
#include "./optimization/multi_extrema.h"

// Complex and quaternion functions
#include "./complex/complex_analysis.h"
#include "./complex/phasor.h"

// Derivative and integral approximation
#include "./calculus/derivation.h"
#include "./calculus/integration.h"
#include "./calculus/odes.h"
#include "./calculus/taylor.h"

// Polynomial class
#include "./polynomial/polynomial.h"
#include "./polynomial/ortho_polyn.h"

// Dual numbers and automatic differentiation
#include "./autodiff/dual.h"
#include "./autodiff/dual_functions.h"
#include "./autodiff/multidual.h"
#include "./autodiff/multidual_functions.h"
#include "./autodiff/dual2.h"
#include "./autodiff/dual2_functions.h"
#include "./autodiff/autodiff.h"

// Pseudorandom number generation
#include "./pseudorandom/pseudorandom.h"
#include "./pseudorandom/prng.h"
#include "./pseudorandom/rand_dist.h"

// Monte Carlo methods
#include "./pseudorandom/montecarlo.h"

#endif

#ifndef THEORETICA_NO_NAMESPACE
/// @namespace th Alias for the theoretica namespace
namespace th = theoretica;
#endif

#endif
