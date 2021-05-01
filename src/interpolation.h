#ifndef UROBORO_INTERP_H
#define UROBORO_INTERP_H

#include "./constants.h"

namespace uroboro {

	// Linear interpolation
	inline real lerp(real number1, real number2, real interp) {
		return (number1 + interp * (number2 - number1));
	}

}

#endif
