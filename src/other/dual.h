#ifndef UROBORO_DUAL_H
#define UROBORO_DUAL_H

#include "../error.h"
#include "../constants.h"

namespace uroboro {


	// Dual number (a + b*epsilon)
	// epsilon such that epsilon^2 = 0
	class dual {
		public:

			real a; // Real part
			real b; // "Dual" part

			// Default constructor, initialize with null values
			dual() : a(0), b(0) {}

			// Initialize from two real numbers
			dual(real real_part, real dual_part)
				: a(real_part), b(dual_part) {}

			~dual() = default;
		
			// Initialize from a vec2
			dual(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
			}

			// Initialize a dual number from a vec2
			inline dual& operator=(const vec2& v) {
				a = v.data[0];
				b = v.data[1];
				return *this;
			}

			// Initialize a dual number from a std::array
			inline dual& operator=(const std::array<real, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			// Return real part
			inline real Re() const {
				return a;
			}

			// Return dual part
			inline real Dual() const {
				return b;
			}

			// Get the inverse of a dual number
			inline dual inverse() const {

				if(a == 0) {
					UMATH_ERROR("dual::inverse", 0, UMATH_ERRCODE::DIV_BY_ZERO);
					return dual(nan(), nan());
				}

				return dual(1.0 / a, -b / square(a));
			}

			// Identity (for consistency)
			inline dual operator+() const {
				return dual(a, b);
			}

			// Sum two dual numbers
			inline dual operator+(const dual& other) const {
				return dual(a + other.a, b + other.b);
			}

			// Sum a real number to a dual number
			inline dual operator+(real r) const {
				return dual(a + r, b);
			}

			// Get the opposite of a dual number
			inline dual operator-() const {
				return dual(-a, -b);
			}

			// Subtract two dual numbers
			inline dual operator-(const dual& other) const {
				return dual(a - other.a, b - other.b);
			}

			// Subtract a real number from a dual number
			inline dual operator-(real r) const {
				return dual(a - r, b);
			}

			// Multiply two dual numbers
			inline dual operator*(const dual& other) const {
				return dual(a * other.a, a * other.b + b * other.a);
			}

			// Multiply a dual number by a real number
			inline dual operator*(real r) const {
				return dual(a * r, b * r);
			}

			// Dual division
			inline dual operator/(const dual& other) const {
				return dual(a / other.a,
					(b * other.a - a * other.b) / square(other.a));
			}

			// Divide a dual number by a real number
			inline dual operator/(real r) const {
				return dual(a / r, b / r);
			}


			// Sum a real number to this dual number
			inline dual& operator+=(real r) {

				a += r;
				return *this;
			}

			// Subtract a real number from this one
			inline dual& operator-=(const dual& other) {

				a -= other.a;
				b -= other.b;
				return *this;
			}

			// Subtract a real number from this dual number
			inline dual& operator-=(real r) {

				a -= r;
				return *this;
			}

			// Multiply this dual number by another one
			inline dual& operator*=(const dual& other) {

				a = (a * other.a);
				b = (a * other.b) + (b * other.a);
				return *this;
			}

			// Multiply this dual number by a real number
			inline dual& operator*=(real r) {

				a *= r;
				b *= r;
				return *this;
			}

			// Divide this dual number by another one
			inline dual& operator/=(const dual& other) {

				a = (a / other.a);
				b = (b * other.a - a * other.b) / square(other.a);
				return *this;
			}

			// Divide a dual number by a real number
			inline dual& operator/=(real r) {

				if(r == 0) {
					UMATH_ERROR("dual::operator/=", 0, UMATH_ERRCODE::DIV_BY_ZERO);
					a = nan();
					b = nan();
					return *this;
				}

				a /= r;
				b /= r;

				return *this;
			}


			// Check whether two dual numbers have the same
			// real and dual parts
			inline bool operator==(const dual& other) {
				return (a == other.a) && (b == other.b);
			}


			// Convert a dual number to a vector
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


			// Convert a dual number to matrix form
			inline mat2 to_mat() const {

				mat2 m;
				m[0][0] = a;
				m[0][1] = 0;
				m[1][0] = a;
				m[1][1] = b;
				return m;
			}

	};


}


#endif
