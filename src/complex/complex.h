
///
/// @file complex.h Complex numbers
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


namespace  theoretica {


	/// @class complex
	/// Complex number in algebraic form \f$a + ib\f$
	template<typename T = real>
	class complex {
		public:

			T a; // Real part
			T b; // Imaginary part

			/// Initialize as \f$(0 + i0)\f$
			complex() : a(0), b(0) {}

			/// Initialize from a real number (zero imaginary part)
			complex(T real_part) : a(real_part), b(0) {}

			/// Initialize from two real numbers
			complex(T real_part, T imag_part) : a(real_part), b(imag_part) {}

			/// Initialize from a real number (zero imaginary part)
			complex(int real_part) : a(real_part), b(0) {}

			/// Initialize a complex number from a std::array
			inline complex& operator=(const std::array<T, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			~complex() = default;

			/// Return real part
			inline T Re() const {
				return a;
			}

			/// Return imaginary part
			inline T Im() const {
				return b;
			}


			/// Extract the real part
			inline friend Re(const complex& z) {
				return z.Re();
			}


			/// Extract the imaginary part
			inline friend T Im(const complex& z) {
				return z.Im();
			}

			/// Get the modulus of a complex number
			inline T modulus() const {
				return sqrt(a * a + b * b);
			}

			/// Get the square modulus of a complex number
			inline T square_modulus() const {
				return a * a + b * b;
			}

			/// Get the complex conjugate of a complex number
			inline complex conjugate() const {
				return complex(a, -b);
			}

			/// Get the inverse of a complex number
			inline complex inverse() const {
				return conjugate() / square_modulus();
			}

			/// Get the argument of a complex number
			inline T arg() const {

				if(b == 0)
					return 0;

				return atan2(b, a);
			}

			/// Identity (for consistency)
			inline complex operator+() const {
				return complex(a, b);
			}

			/// Sum two complex numbers
			inline complex operator+(const complex& other) const {
				return complex(a + other.a, b + other.b);
			}

			/// Sum a T number to a complex number
			inline complex operator+(T r) const {
				return complex(a + r, b);
			}

			/// Get the opposite of a complex number
			inline complex operator-() const {
				return complex(-a, -b);
			}

			/// Subtract two complex numbers
			inline complex operator-(const complex& other) const {
				return complex(a - other.a, b - other.b);
			}

			/// Subtract a real number from a complex number
			inline complex operator-(T r) const {
				return complex(a - r, b);
			}

			/// Multiply two complex numbers
			inline complex operator*(const complex& other) const {
				return complex((a * other.a) - (b * other.b), (a * other.b) + (b * other.a));
			}

			/// Multiply a complex number by a real number
			inline complex operator*(T r) const {
				return complex(a * r, b * r);
			}

			/// Complex division
			inline complex operator/(const complex& other) const {
				T m = other.square_modulus();
				return complex((a * other.a + b * other.b) / m,
								(b * other.a - a * other.b) / m);
			}

			/// Divide a complex number by a real number
			inline complex operator/(T r) const {
				return complex(a / r, b / r);
			}


			/// Add a real number to this one
			inline complex& operator+=(const complex& other) {

				a += other.a;
				b += other.b;
				return *this;
			}


			/// Sum a real number to this complex number
			inline complex& operator+=(T r) {

				a += r;
				return *this;
			}

			/// Subtract a real number from this one
			inline complex& operator-=(const complex& other) {

				a -= other.a;
				b -= other.b;
				return *this;
			}

			/// Subtract a real number from this complex number
			inline complex& operator-=(T r) {

				a -= r;
				return *this;
			}

			/// Multiply this complex number by another one
			inline complex& operator*=(const complex& other) {

				a = (a * other.a) - (b * other.b);
				b = (a * other.b) + (b * other.a);
				return *this;
			}

			/// Multiply this complex number by a real number
			inline complex& operator*=(T r) {

				a *= r;
				b *= r;
				return *this;
			}

			/// Divide this complex number by another one
			inline complex& operator/=(const complex& other) {

				T m = other.square_modulus();
				a = (a * other.a + b * other.b) / m;
				b = (b * other.a - a * other.b) / m;
				return *this;
			}

			/// Divide a complex number by a real number
			inline complex& operator/=(T r) {

				if(r == 0) {
					TH_MATH_ERROR("complex::operator/=", r, DIV_BY_ZERO);
					a = nan();
					b = nan();
					return *this;
				}

				a /= r;
				b /= r;

				return *this;
			}


			/// Check whether two complex numbers have the same
			/// real and imaginary parts
			inline bool operator==(const complex& other) {
				return (a == other.a) && (b == other.b);
			}


			/// Check whether two complex numbers have the same
			/// real and imaginary parts
			inline bool operator!=(const complex& other) {
				return !(*this == other);
			}


			/// Convert a complex number to a vector
			template<typename Vector>
			inline Vector to_vec() const {
				Vector res;
				res[0] = a;
				res[1] = b;
				return res;
			}

			/// Initialize from a vector
			template<typename Vector>
			inline void from_vec(const Vector& v) {
				a = v[0];
				b = v[1];
			}


			/// Construct a complex number representing a rotation
			/// of <rad> radians in 2 dimensions
			inline static complex rotor(T rad) {
				return complex(cos(rad), sin(rad));
			}


			inline static complex i() {
				return complex(0, 1);
			}


			// Friend operators to enable equations of the form
			// (real) op. (complex)

			inline friend complex operator+(T r, const complex& z) {
				return z + r;
			}

			inline friend complex operator-(T r, const complex& z) {
				return -z + r;
			}

			inline friend complex operator*(T r, const complex& z) {
				return z * r;
			}

			inline friend complex operator/(T r, const complex& z) {
				return complex(r, 0) / z;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the complex number to string representation
			inline std::string to_string() const {

				std::stringstream res;

				res << a;
				res << (b >= 0 ? " + " : " - ");
				
				if(abs(b) != 1)
					res << abs(b);

				res << "i";

				return res.str();
			}


			/// Stream the complex number in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const complex& obj) {
				return out << obj.to_string();
			}

#endif

	};


	// To-DO: Does not compile on GCC v6	
	// constexpr complex operator""i(long double b) {
	// 	return complex(0, b);
	// }

}

#endif
