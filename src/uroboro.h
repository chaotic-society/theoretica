#ifndef UROBORO_H
#define UROBORO_H

#include "./constants.h"
#include "./error.h"
#include "./real_analysis.h"
#include "./algebra/vec.h"
#include "./algebra/mat.h"
#include "./complex/complex.h"
#include "./complex/quat.h"

#ifdef UROBORO_INCLUDE_ALL
#include "./vec_buff.h"
#include "./utility.h"
#include "./other/interpolation.h"
#include "./statistics/statistics.h"
#include "./approx.h"
#include "./complex/complex_analysis.h"
#include "./complex/phasor.h"
#include "./calculus/derivation.h"
#include "./calculus/integration.h"
#include "./polynomial.h"
#include "./statistics/distributions.h"
#include "./autodiff/dual.h"
#include "./autodiff/dual_functions.h"
#endif

#ifndef UROBORO_NO_NAMESPACE
namespace umath = uroboro;
#endif

#endif
