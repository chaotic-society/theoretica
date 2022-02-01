#ifndef UROBORO_PHASOR
#define UROBORO_PHASOR

#include "../real_analysis.h"

namespace uroboro {

	class phasor {
		public:

			real phase;
			real modulus;

			// Initialize as 0/0
			phasor() : phase(0), modulus(0) {}

			// Initialize from two real numbers
			phasor(real modulus, real phase) {
				this->phase = phase;
				this->modulus = modulus;
			}

			// Initialize from another phasor
			phasor(const phasor& other) {
				phase = other.phase;
				modulus = other.modulus;
			}

			// Initialize from a complex number
			phasor(const complex& z) {
				phase = z.arg();
				modulus = z.modulus();
			}

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


			// TO-DO +=, *= ... operators


			// Transform a phasor to a complex number
			complex to_complex() const {
				return complex(
					modulus * cos(phase),
					modulus * sin(phase));
			}

	};

}

#endif
