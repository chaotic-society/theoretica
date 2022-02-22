#ifndef UROBORO_PHASOR
#define UROBORO_PHASOR

#include "../real_analysis.h"
#include "./complex.h"


namespace uroboro {

	// Complex number in exponential (modulus * e^(i * phase))
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

				if(other.modulus == 0) {
					UMATH_ERROR("phasor::operator/", other.modulus, DIV_BY_ZERO);
					return phasor(nan(), nan());
				}
					

				return phasor(
					modulus / other.modulus,
					phase - other.phase);
			}


			// Add a phasor to this one
			inline phasor& operator+=(const phasor& other) {

				if(phase == other.phase) {
					modulus += other.modulus;
				} else {
					*this = phasor(to_complex() + other.to_complex());
				}

				return *this;
			}


			// Subtract a phasor from this one
			inline phasor& operator-=(const phasor& other) {

				if(phase == other.phase)
					modulus -= other.modulus;
				else
					*this = phasor(to_complex() - other.to_complex());

				return *this;
			}


			// Multiply this phasor by another one
			inline phasor& operator*=(const phasor& other) {

				modulus *= other.modulus;
				phase += other.phase;

				return *this;
			}


			// Divide this phasor by another one
			inline phasor& operator/=(const phasor& other) {

				if(other.modulus == 0) {
					UMATH_ERROR("phasor::operator/=", other.modulus, DIV_BY_ZERO);
					modulus = nan();
					phase = nan();
					return *this;
				}

				modulus /= other.modulus;
				phase -= other.phase;

				return *this;
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
