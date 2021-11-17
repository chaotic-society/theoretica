#ifndef UROBORO_H
#define UROBORO_H

#include "./common.h"
#include "./vec.h"
#include "./mat.h"
#include "./complex/complex.h"
#include "./complex/quat.h"
#include "./complex/phasor.h"

#ifdef UROBORO_INCLUDE_ALL
#include "./vec_buff.h"
#include "./utility.h"
#include "./interpolation.h"
#include "./statistics/statistics.h"
#include "./approx.h"
#include "./complex/complex_functions.h"
#include "./derivation.h"
#include "./integration.h"
#include "./polynomial.h"
#include "./statistics/distributions.h"
#endif

#ifndef UROBORO_NO_NAMESPACE
namespace umath = uroboro;
#endif

#endif
