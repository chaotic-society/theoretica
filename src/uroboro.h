#ifndef UROBORO_H
#define UROBORO_H

#include "./common.h"
#include "./vec.h"
#include "./mat.h"
#include "./complex.h"
#include "./quat.h"
#include "./phasor.h"

#ifdef UROBORO_INCLUDE_ALL
#include "./vec_buff.h"
#include "./utility.h"
#include "./interpolation.h"
#include "./statistics.h"
#include "./approx.h"
#include "./complex_functions.h"
#include "./derivation.h"
#include "./integration.h"
#include "./polynomial.h"
#include "./distributions.h"
#endif

#ifndef UROBORO_NO_NAMESPACE
namespace umath = uroboro;
#endif

#endif
