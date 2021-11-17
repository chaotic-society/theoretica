#ifndef UROBORO_PHASOR
#define UROBORO_PHASOR

#include "../common.h"

namespace uroboro {

	class phasor {
		public:

			real phase;
			real modulus;

			phasor(real modulus, real phase) {
				this->phase = phase;
				this->modulus = modulus;
			}


			phasor(const phasor& other) {
				phase = other.phase;
				modulus = other.modulus;
			}


			phasor(const complex& z) {
				phase = z.arg();
				modulus = z.modulus();
			}


			phasor() = default;
			~phasor() = default;


			// Add two phasors
			inline phasor operator+(const phasor& other) const {

				if(phase == other.phase)
					return phasor(modulus + other.modulus, phase);

				return phasor(to_complex() + other.to_complex());
			}


			// Subtract two phasors
			inline phasor operator-(const phasor& other) const {

				if(phase == other.phase)
					return phasor(modulus - other.modulus, phase);

				return phasor(to_complex() - other.to_complex());
			}


			// Multiply two phasors
			inline phasor operator*(const phasor& other) const {
				return phasor(
					modulus * other.modulus,
					phase + other.phase);
			}


			// Divide two phasors
			inline phasor operator/(const phasor& other) const {

				if(other.modulus == 0)
					// throw ...
					return phasor(0, 0);

				return phasor(
					modulus / other.modulus,
					phase - other.phase);
			}


			// Transform a phasor to a complex number
			complex to_complex() const {
				return complex(
					modulus * cos(phase),
					modulus * sin(phase));
			}

	};

}

#endif
