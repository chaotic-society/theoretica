#ifndef UROBORO_COMPLEX_H
#define UROBORO_COMPLEX_H

#include "../real_analysis.h"
#include "../algebra/vec.h"
#include <array>

namespace  uroboro {

	// Complex number in algebraic form (a + ib)
	class complex {
		public:

			real a; // Real part
			real b; // Imaginary part

			// Initialize as (0 + i0)
			complex() : a(0), b(0) {}

			// Initialize from two real numbers
			complex(real r, real i) : a(r), b(i) {}

			// Initialize from a vec2
			complex(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
			}

			// Initialize a complex number from a vec2
			inline complex& operator=(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
				return *this;
			}

			// Initialize a complex number from a std::array
			inline complex& operator=(const std::array<real, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			~complex() = default;

			// Return real part
			inline real Re() const {
				return a;
			}

			// Return imaginary part
			inline real Im() const {
				return b;
			}

			// Get the modulus of a complex number
			inline real modulus() const {
				return uroboro::sqrt(a * a + b * b);
			}

			// Get the square modulus of a complex number
			inline real square_modulus() const {
				return a * a + b * b;
			}

			// Get the complex conjugate of a complex number
			inline complex conjugate() const {
				return complex(a, -b);
			}

			// Get the inverse of a complex number
			inline complex inverse() const {
				return conjugate() / square_modulus();
			}

			// Get the argument of a complex number
			inline real arg() const {

				if(b == 0)
					return 0;

				return atan2(b, a);
			}

			// Identity (for consistency)
			inline complex operator+() const {
				return complex(a, b);
			}

			// Sum two complex numbers
			inline complex operator+(const complex& other) const {
				return complex(a + other.a, b + other.b);
			}

			// Sum a real number to a complex number
			inline complex operator+(real r) const {
				return complex(a + r, b);
			}

			// Get the opposite of a complex number
			inline complex operator-() const {
				return complex(-a, -b);
			}

			// Subtract two complex numbers
			inline complex operator-(const complex& other) const {
				return complex(a - other.a, b - other.b);
			}

			// Subtract a real number from a complex number
			inline complex operator-(real r) const {
				return complex(a - r, b);
			}

			// Multiply two complex numbers
			inline complex operator*(const complex& other) const {
				return complex((a * other.a) - (b * other.b), (a * other.b) + (b * other.a));
			}

			// Multiply a complex number by a real number
			inline complex operator*(real r) const {
				return complex(a * r, b * r);
			}

			// Complex division
			inline complex operator/(const complex& other) const {
				real m = other.square_modulus();
				return complex((a * other.a + b * other.b) / m,
								(b * other.a - a * other.b) / m);
			}

			// Divide a complex number by a real number
			inline complex operator/(real r) const {
				return complex(a / r, b / r);
			}


			// Sum a real number to this complex number
			inline complex& operator+=(real r) {

				a += r;
				return *this;
			}

			// Subtract a real number from this one
			inline complex& operator-=(const complex& other) {

				a -= other.a;
				b -= other.b;
				return *this;
			}

			// Subtract a real number from this complex number
			inline complex& operator-=(real r) {

				a -= r;
				return *this;
			}

			// Multiply this complex number by another one
			inline complex& operator*=(const complex& other) {

				a = (a * other.a) - (b * other.b);
				b = (a * other.b) + (b * other.a);
				return *this;
			}

			// Multiply this complex number by a real number
			inline complex& operator*=(real r) {

				a *= r;
				b *= r;
				return *this;
			}

			// Divide this complex number by another one
			inline complex& operator/=(const complex& other) {

				real m = other.square_modulus();
				a = (a * other.a + b * other.b) / m;
				b = (b * other.a - a * other.b) / m;
				return *this;
			}

			// Divide a complex number by a real number
			inline complex& operator/=(real r) {

				if(r == 0) {
					// throw ...
					a = 0;
					b = 0;
				}

				a /= r;
				b /= r;

				return *this;
			}


			// Check whether two complex numbers have the same
			// real and imaginary parts
			inline bool operator==(const complex& other) {
				return (a == other.a) && (b == other.b);
			}


			// Convert a complex number to a vector
			inline vec2 to_vec() const {
				vec2 res;
				res.data[0] = a;
				res.data[1] = b;
				return res;
			}

			// Initialize from a vector
			inline void from_vec(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
			}


			// Construct a complex number representing a rotation
			// of <rad> radians in 2 dimensions
			inline static complex rotor(real rad) {
				return complex(cos(rad), sin(rad));
			}


			inline static complex i() {
				return complex(0, 1);
			}

	};

}

#endif
