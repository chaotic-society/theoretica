
///
/// @file phasor.h Phasor (complex number in exponential form)
///

#ifndef THEORETICA_PHASOR
#define THEORETICA_PHASOR

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/real_analysis.h"
#include "./complex.h"


namespace theoretica {


	/// @class phasor
	/// Complex number in exponential form \f$\rho e^{i \theta}\f$
	template<typename Type = real>
	class phasor {
		public:

			/// Modulus of the complex number
			Type modulus;

			/// Phase of the complex number
			Type phase;


			/// Initialize as 0/0
			phasor() : phase(0), modulus(0) {}


			/// Initialize from two real numbers
			phasor(Type modulus, Type phase) {
				this->phase = phase;
				this->modulus = modulus;
			}


			/// Initialize from another phasor
			phasor(const phasor& other) {
				phase = other.phase;
				modulus = other.modulus;
			}

			/// Construct a phasor from a complex number
			template<typename T>
			phasor(const complex<T>& z) {
				phase = z.arg();
				modulus = z.norm();
			}

			/// Construct a phasor from a real number
			phasor(Type r) {
				modulus = abs(r);
				phase = (r >= 0 ? 0 : PI);
			}


			/// Assignment operator
			inline phasor& operator=(const phasor& z) {
				phase = z.phase;
				modulus = z.modulus;
				return *this;
			}


			/// Assignment operator from a 2D array
			/// as {modulus, phase}.
			template<typename T>
			inline phasor& operator=(const std::array<T, 2>& v) {
				modulus = v[0];
				phase = v[1];
				return *this;
			}


			~phasor() = default;


			/// Get the real part of the complex number
			inline Type Re() const {
				return modulus * cos(phase);
			}


			/// Extract the real part of the complex number
			inline friend Type Re(const phasor& z) {
				return z.Re();
			}


			/// Get the imaginary part of the complex number
			inline Type Im() const {
				return modulus * sin(phase);
			}


			/// Extract the imaginary part of the complex number
			inline friend Type Im(const phasor& z) {
				return z.Im();
			}


			/// Compute the conjugate of the complex number
			inline phasor conjugate() const {
				return phasor(modulus, -phase);
			}


			/// Compute the square norm of the complex number
			inline Type sqr_norm() const {
				return modulus * modulus;
			}


			/// Compute the norm of the complex number
			inline Type norm() const {
				return modulus;
			}


			/// Compute the inverse of the complex number
			inline Type inverse() const {

				if(modulus < MACH_EPSILON) {
					TH_MATH_ERROR("phasor::inverse", modulus, DIV_BY_ZERO);
					return (Type) nan();
				}

				return phasor(
					((Type) 1.0) / modulus,
					-phase
				);
			}


			/// Invert the complex number
			inline phasor& invert() {

				if(modulus < MACH_EPSILON) {
					TH_MATH_ERROR("phasor::invert", modulus, DIV_BY_ZERO);
					modulus = (Type) nan();
					phase = (Type) nan();
					return *this;
				}

				phase = phase * -1;
				modulus = ((Type) 1.0) / modulus;
				return *this;
			}


			/// Get the argument of the complex number
			inline Type arg() const {
				return phase;
			}


			/// Add two phasors
			/// @note This operation is particularly slow
			/// for phasors as opposed to complex numbers
			/// in algebraic form.
			inline phasor operator+(const phasor& other) const {

				if(abs(phase - other.phase) < MACH_EPSILON)
					return phasor(modulus + other.modulus, phase);

				return phasor(to_complex() + other.to_complex());
			}


			/// Subtract two phasors
			/// @note This operation is particularly slow
			/// for phasors as opposed to complex numbers
			/// in algebraic form.
			inline phasor operator-(const phasor& other) const {

				if(abs(phase - other.phase) < MACH_EPSILON)
					return phasor(modulus - other.modulus, phase);

				return phasor(to_complex() - other.to_complex());
			}


			/// Multiply two phasors
			inline phasor operator*(const phasor& other) const {
				return phasor(
					modulus * other.modulus,
					phase + other.phase);
			}


			/// Multiply a complex number in algebraic form
			/// and a phasor
			template<typename T>
			inline phasor operator*(const complex<T>& other) const {
				return phasor(
					modulus * other.norm(),
					phase + other.arg());
			}


			/// Multiply a complex number in algebraic form
			/// and a phasor
			template<typename T>
			inline friend phasor operator*(const complex<T>& z, const phasor& w) {
				return phasor(
					z.norm() * w.modulus,
					z.arg() + w.phase
				);
			}


			/// Divide two phasors
			inline phasor operator/(const phasor& other) const {

				if(abs(other.modulus) < MACH_EPSILON) {
					TH_MATH_ERROR("phasor::operator/", other.modulus, DIV_BY_ZERO);
					return phasor(nan(), nan());
				}
					

				return phasor(
					modulus / other.modulus,
					phase - other.phase);
			}


			/// Add a phasor to this one
			inline phasor& operator+=(const phasor& other) {

				if(phase == other.phase) {
					modulus += other.modulus;
				} else {
					*this = phasor(to_complex() + other.to_complex());
				}

				return *this;
			}


			/// Subtract a phasor from this one
			inline phasor& operator-=(const phasor& other) {

				if(phase == other.phase)
					modulus -= other.modulus;
				else
					*this = phasor(to_complex() - other.to_complex());

				return *this;
			}


			/// Multiply this phasor by another one
			inline phasor& operator*=(const phasor& other) {

				modulus *= other.modulus;
				phase += other.phase;

				return *this;
			}


			/// Divide this phasor by another one
			inline phasor& operator/=(const phasor& other) {

				if(abs(other.modulus) < MACH_EPSILON) {
					TH_MATH_ERROR("phasor::operator/=", other.modulus, DIV_BY_ZERO);
					modulus = nan();
					phase = nan();
					return *this;
				}

				modulus /= other.modulus;
				phase -= other.phase;

				return *this;
			}


			/// Check whether two phasors are the same
			inline bool operator==(const phasor& z) const {
				return (modulus == z.modulus) && (phase == z.phase);
			}


			/// Check whether two phasors are not the same
			inline bool operator!=(const phasor& z) const {
				return !(*this == z);
			}


			/// Transform a phasor to a complex number
			template<typename T = Type>
			complex<T> to_complex() const {
				return complex<T>(
					modulus * cos(phase),
					modulus * sin(phase));
			}


			/// Cast to complex
			template<typename T>
			inline operator complex<T> () {
				return to_complex<T>();
			}


			/// Construct a phasor representing a rotation
			/// of <rad> radians in 2 dimensions
			inline static phasor rotor(Type rad) {
				return phasor((Type) 1.0, rad);
			}


			/// Imaginary unit in exponential form
			inline static phasor i() {
				return phasor(1, PI / 2.0);
			}


			// Friend operators to enable equations of the form
			// (real) op. (phasor)

			inline friend phasor operator+(Type r, const phasor& z) {
				return z + phasor(r);
			}

			inline friend phasor operator-(Type r, const phasor& z) {
				return phasor(z.modulus, z.phase + PI) + phasor(r);
			}

			inline friend phasor operator*(Type r, const phasor& z) {
				return z * phasor(r);
			}

			inline friend phasor operator/(Type r, const phasor& z) {
				return phasor(r) / z;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the phasor to string representation
			inline std::string to_string(const std::string& separator = ", ") const {

				std::stringstream res;
				res << modulus << "/" << phase;
				return res.str();
			}


			/// Stream the phasor in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const phasor& obj) {
				return out << obj.to_string();
			}

#endif

	};

}

#endif
