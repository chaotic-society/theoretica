#ifndef UROBORO_INTERP_H
#define UROBORO_INTERP_H

#include "common.h"

namespace uroboro {

	inline real lerp(real number1, real number2, real interp) {
		return (number1 + interp * (number2 - number1));
	}

	inline real smoothstep(real number1, real number2, real interp) {
		const real result = clamp((interp - number1) / (number2 - number1), 0.0, 1.0);
		return result * result * (3 - 2 * result);
	}

	inline real smootherstep(real number1, real number2, real interp) {
		const real result = clamp((interp - number1) / (number2 - number1), 0.0, 1.0);
		return result * result * result * (result * (result * 6 - 15) + 10);
	}


}

#endif
