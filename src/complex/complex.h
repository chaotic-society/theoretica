
///
/// @file complex.h Complex number class
///

#ifndef THEORETICA_COMPLEX_H
#define THEORETICA_COMPLEX_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include <array>

#include "../core/error.h"
#include "../core/real_analysis.h"


namespace theoretica {


	/// @class complex
	/// Complex number in algebraic form \f$a + ib\f$
	template<typename Type = real>
	class complex {
		public:

			/// Real part
			Type a;

			/// Imaginary part
			Type b;


			/// Default constructor, initializes
			/// the number to zero
			complex() : a(0), b(0) {}


			/// Construct a complex number from its real
			/// and complex parts.
			complex(Type real_part, Type imag_part)
				: a(real_part), b(imag_part) {}


			/// Construct a complex number from a real number,
			/// with zero imaginary part.
			complex(Type real_part) : a(real_part), b(0) {}


			/// Copy constructor
			complex(const complex& z) : a(z.a), b(z.b) {}


			/// Assignment operator
			inline complex& operator=(const complex& z) {
				a = z.a;
				b = z.b;
				return *this;
			}


			/// Assignment operator from a 2D array
			template<typename T = Type>
			inline complex& operator=(const std::array<T, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}


			/// Get the real part of the complex number
			inline Type Re() const {
				return a;
			}


			/// Get the real part of the complex number
			inline Type& Re() {
				return a;
			}


			/// Extract the real part of the complex number
			inline friend Type Re(const complex& z) {
				return z.a;
			}


			/// Extract the real part of the complex number
			inline friend Type& Re(complex& z) {
				return z.a;
			}


			/// Get the imaginary part of the complex number
			inline Type Im() const {
				return b;
			}


			/// Get the imaginary part of the complex number
			inline Type& Im() {
				return b;
			}


			/// Extract the imaginary part of the complex number
			inline friend Type Im(const complex& z) {
				return z.b;
			}


			/// Extract the imaginary part of the complex number
			inline friend Type& Im(complex& z) {
				return z.b;
			}


			/// Compute the conjugate of the complex number
			inline complex conjugate() const {
				return complex(a, -b);
			}


			/// Compute the square norm of the complex number
			inline Type sqr_norm() const {
				return a * a + b * b;
			}


			/// Compute the norm of the complex number
			inline Type norm() const {
				return sqrt(sqr_norm());
			}


			/// Compute the inverse of the complex number
			inline complex inverse() const {

				const Type n = sqr_norm();

				if(n < MACH_EPSILON) {
					TH_MATH_ERROR("complex::inverse", n, MathError::DivByZero);
					return complex(nan(), nan());
				}

				return conjugate() / n;
			}


			/// Invert the complex number
			inline complex& invert() {

				const Type n = sqr_norm();

				if(n < MACH_EPSILON) {
					TH_MATH_ERROR("complex::invert", n, MathError::DivByZero);
					a = (Type) nan();
					b = (Type) nan();
					return *this;
				}

				a = a / n;
				b = -b / n;

				return *this;
			}


			/// Get the argument of the complex number
			inline Type arg() const {

				// Real number case
				if(abs(b) < MACH_EPSILON) {
					return (a >= 0) ? ((Type) 0) : ((Type) PI);
				}

				// Pure imaginary number case
				if(abs(a) < MACH_EPSILON) {
					return (b >= 0) ? ((Type) PI / 2.0) : ((Type) -PI / 2.0);
				}

				// Use the 2-parameter arctangent in the general case
				return atan2(b, a);
			}


			/// Identity (for consistency)
			inline complex operator+() const {
				return complex(a, b);
			}


			/// Add two complex numbers
			inline complex operator+(const complex& z) const {
				return complex(a + z.a, b + z.b);
			}


			/// Get the opposite of the complex number
			inline complex operator-() const {
				return complex(-a, -b);
			}


			/// Subtract two complex numbers
			inline complex operator-(const complex& z) const {
				return complex(a - z.a, b - z.b);
			}


			/// Multiply two complex numbers
			inline complex operator*(const complex& z) const {
				return complex(
					a * z.a - b * z.b,
					a * z.b + b * z.a
				);
			}


			/// Divide two complex numbers
			inline complex operator/(const complex& z) const {
				return operator*(z.inverse());
			}


			/// Add a complex number to this one
			inline complex& operator+=(const complex& z) {
				return (*this = complex(a + z.a, b + z.b));
			}


			/// Subtract a complex number from this one
			inline complex& operator-=(const complex& z) {
				return (*this = complex(a - z.a, b - z.b));
			}


			/// Multiply the complex number by another
			inline complex& operator*=(const complex& z) {
				return (*this = complex(
					a * z.a - b * z.b,
					a * z.b + b * z.a
				));
			}


			/// Divide the complex number by another
			inline complex& operator/=(const complex& z) {
				return (*this = operator*(z.inverse()));
			}


			/// Add a real number to the complex number
			inline complex operator+(Type k) const {
				return complex(a + k, b);
			}


			/// Subtract a real number from the complex number
			inline complex operator-(Type k) const {
				return complex(a - k, b);
			}


			/// Multiply the complex number by a real number
			inline complex operator*(Type k) const {
				return complex(a * k, b * k);
			}


			/// Divide the complex number by a real number
			inline complex operator/(Type k) const {

				if(abs(k) < MACH_EPSILON) {
					TH_MATH_ERROR("complex::operator/", k, MathError::DivByZero);
					return complex(nan(), nan());
				}

				return complex(a / k, b / k);
			}


			/// Add a real number to the complex number
			inline complex& operator+=(Type k) {
				return (*this = complex(a + k, b));
			}


			/// Subtract a real number from the complex number
			inline complex& operator-=(Type k) {
				return (*this = complex(a - k, b));
			}


			/// Multiply the complex number by a real number
			inline complex& operator*=(Type k) {
				return (*this = complex(a * k, b * k));
			}


			/// Divide the complex number by a real number
			inline complex& operator/=(Type k) {

				if(abs(k) < MACH_EPSILON) {
					TH_MATH_ERROR("complex::operator/=", k, MathError::DivByZero);
					return complex(nan(), nan());
				}

				return (*this = complex(a / k, b / k));
			}


			/// Check whether two complex numbers are the same
			inline bool operator==(const complex& z) const {
				return (a == z.a) && (b == z.b);
			}


			/// Check whether two complex numbers are not the same
			inline bool operator!=(const complex& z) const {
				return !(*this == z);
			}


			/// Construct a complex number representing a rotation
			/// of <rad> radians in 2 dimensions
			inline static complex rotor(Type rad) {
				return complex(cos(rad), sin(rad));
			}


			/// Imaginary unit
			inline static complex i() {
				return complex(0, 1);
			}


			// Friend operators to enable equations of the form
			// (real) op. (complex)

			inline friend complex operator+(Type r, const complex& z) {
				return z + r;
			}

			inline friend complex operator-(Type r, const complex& z) {
				return -z + r;
			}

			inline friend complex operator*(Type r, const complex& z) {
				return z * r;
			}

			inline friend complex operator/(Type r, const complex& z) {
				return complex(r, 0) / z;
			}


			/// Narrowing cast to a real number. If the imaginary part
			/// is greater in module than the machine epsilon, NaN is returned.
			inline operator real () {

				if (abs(b) >= MACH_EPSILON) {
					return nan();
				}

				return a;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the complex number to string representation
			inline std::string to_string() const {

				std::stringstream res;

				// TO-DO Check whether Type has comparison
				// operators (bicomplex with Type=complex !)

				res << a;
				res << (b >= 0 ? " + " : " - ");
				
				if(abs(b) != 1)
					res << abs(b);

				res << "i";

				return res.str();
			}


			/// Convert the complex number to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the complex number in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const complex& obj) {
				return out << obj.to_string();
			}

#endif

	};

}

#endif
