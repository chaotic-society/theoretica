
///
/// @file uroboro.h General include file
/// This file will include core header files.
/// If `UROBORO_INCLUDE_ALL` is defined, it will include all header files of the library
///

#ifndef UROBORO_H
#define UROBORO_H

// Core headers
#include "./constants.h"
#include "./error.h"
#include "./real_analysis.h"
#include "./algebra/vec.h"
#include "./algebra/mat.h"
#include "./complex/complex.h"
#include "./complex/quat.h"

// All other headers
#ifdef UROBORO_INCLUDE_ALL
#include "./vec_buff.h"
#include "./utility.h"
#include "./interpolation/spline_interp.h"
#include "./interpolation/polyn_interp.h"
#include "./statistics/statistics.h"
#include "./approx/roots.h"
#include "./approx/extrema.h"
#include "./complex/complex_analysis.h"
#include "./complex/phasor.h"
#include "./calculus/derivation.h"
#include "./calculus/integration.h"
#include "./polynomial/polynomial.h"
#include "./statistics/distributions.h"
#include "./autodiff/dual.h"
#include "./autodiff/dual_functions.h"
#include "./autodiff/autodiff.h"
#include "./pseudo_random/pseudo_random_algo.h"
#endif

/// @namespace umath Alias for the uroboro namespace
#ifndef UROBORO_NO_NAMESPACE
namespace umath = uroboro;
#endif

#endif
