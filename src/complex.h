#ifndef UROBORO_COMPLEX_H
#define UROBORO_COMPLEX_H

#include "./common.h"
#include "./vec.h"
#include <array>

namespace  uroboro {

	class complex {
		public:

			real a; // Real part
			real b; // Imaginary part

			complex() {
				a = 0;
				b = 0;
			}

			complex(real r, real i) {
				a = r;
				b = i;
			}

			complex(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
			}

			inline complex& operator=(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
				return *this;
			}

			inline complex& operator=(const std::array<real, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			~complex() {}

			inline real Re() const {
				return a;
			}

			inline real Im() const {
				return b;
			}

			inline real modulus() const {
				return uroboro::sqrt(a * a + b * b);
			}

			inline real square_modulus() const {
				return a * a + b * b;
			}

			inline real magnitude() const {
				return modulus();
			}

			inline real square_magnitude() const {
				return square_modulus();
			}

			inline complex conjugate() const {
				return complex(a, -b);
			}

			inline complex inverse() const {
				return conjugate() / square_modulus();
			}

			inline complex operator+(const complex& other) const {
				return complex(a + other.a, b + other.b);
			}

			inline complex operator-(const complex& other) const {
				return complex(a - other.a, b - other.b);
			}

			inline complex operator*(const complex& other) const {
				return complex((a * other.a) - (b * other.b), (a * other.b) + (b * other.a));
			}

			inline complex operator/(real r) const {
				return complex(a / r, b / r);
			}

			inline complex operator/(const complex& other) const {
				real m = other.square_modulus();
				return complex((a * other.a + b * other.b) / m,
								(b * other.a - a * other.b) / m);
			}

			// operator += -= *= /= ...

			inline bool operator==(const complex& other) {
				return (a == other.a) && (b == other.b);
			}

			inline vec2 to_vec() const {
				vec2 res;
				res.data[0] = a;
				res.data[1] = b;
				return res;
			}

			inline void from_vec(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
			}


			inline static complex rotor(real rad) {
				return complex(uroboro::cos(rad), uroboro::sin(rad));
			}

	};

}

#endif
