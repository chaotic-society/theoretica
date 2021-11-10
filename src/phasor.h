#ifndef UROBORO_PHASOR
#define UROBORO_PHASOR

#include "./common.h"

namespace uroboro {

	class phasor {
		public:

			real phase;
			real modulus;

			phasor(real phase, real modulus) {
				this->phase = phase;
				this->modulus = modulus;
			}


			phasor() = default;
			~phasor() = default;

			complex to_complex() {
				return complex(
					modulus * sin(phase),
					modulus * cos(phase));
			}

	};

}

#endif
