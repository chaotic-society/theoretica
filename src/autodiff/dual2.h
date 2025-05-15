
///
/// @file dual2.h Second order dual number class
///

#ifndef THEORETICA_DUAL2_H
#define THEORETICA_DUAL2_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/error.h"
#include "../core/constants.h"
#include "../algebra/algebra_types.h"


namespace theoretica {


	///
	/// @class dual2
	/// Second order dual number class.
	/// Implemented as \f$a + b * \epsilon_1 + c * \epsilon_2\f$
	/// where epsilon is such that \f$\epsilon_1^3 = 0\f$
	/// and \f$\epsilon_1 * \epsilon_2 = 0\f$
	///
	class dual2 {
		public:

			real a; // Real part
			real b; // First order dual part
			real c; // Second order dual part

			/// Default constructor, initialize with null values
			dual2() : a(0), b(0), c(0) {}

			/// Initialize from two real numbers
			dual2(real real_part, real dual1_part, real dual2_part)
				: a(real_part), b(dual1_part), c(dual2_part) {}

			/// Initialize from two real numbers
			dual2(real real_part, real dual1_part)
				: a(real_part), b(dual1_part), c(0) {}

			/// Initialize from a real number
			dual2(real real_part)
				: a(real_part), b(0), c(0) {}

			~dual2() = default;
		
			/// Initialize from a vec3
			dual2(const vec3& v) {
				a = v[0];
				b = v[1];
				c = v[2];
			}

			/// Initialize a dual number from a vec3
			inline dual2& operator=(const vec3& v) {
				a = v[0];
				b = v[1];
				c = v[2];
				return *this;
			}

			/// Initialize a dual number from a real number
			inline dual2& operator=(real x) {
				a = x;
				b = 0;
				c = 0;
				return *this;
			}

			/// Initialize a dual number from a std::array
			inline dual2& operator=(const std::array<real, 3>& v) {
				a = v[0];
				b = v[1];
				c = v[2];
				return *this;
			}

			/// Return real part
			inline real Re() const {
				return a;
			}

			/// Return first order dual part
			inline real Dual1() const {
				return b;
			}

			/// Return second order dual part
			inline real Dual2() const {
				return c;
			}

			/// Get the dual conjugate
			inline dual2 conjugate() const {
				return dual2(a, -b, -c);
			}

			/// Get the inverse of a dual number
			inline dual2 inverse() const {

				if(a == 0) {
					TH_MATH_ERROR("dual2::inverse", 0, DIV_BY_ZERO);
					return dual2(nan(), nan(), nan());
				}

				return dual2(1.0 / a, -b / square(a), 2 * b / cube(a));
			}

			/// Identity (for consistency)
			inline dual2 operator+() const {
				return dual2(a, b, c);
			}

			/// Sum two dual numbers
			inline dual2 operator+(const dual2& other) const {
				return dual2(a + other.a, b + other.b, c + other.c);
			}

			/// Sum a real number to a dual number
			inline dual2 operator+(real r) const {
				return dual2(a + r, b, c);
			}

			/// Get the opposite of a dual number
			inline dual2 operator-() const {
				return dual2(-a, -b, -c);
			}

			/// Subtract two dual numbers
			inline dual2 operator-(const dual2& other) const {
				return dual2(a - other.a, b - other.b, c - other.c);
			}

			/// Subtract a real number from a dual number
			inline dual2 operator-(real r) const {
				return dual2(a - r, b, c);
			}

			/// Multiply two dual numbers
			inline dual2 operator*(const dual2& other) const {
				return dual2(a * other.a,
					a * other.b + b * other.a,
					a * other.c + 2 * b * other.b + c * other.a);
			}

			/// Multiply a dual number by a real number
			inline dual2 operator*(real r) const {
				return dual2(a * r, b * r, c * r);
			}

			/// Dual division
			inline dual2 operator/(const dual2& other) const {
				return operator*(other.inverse());
			}

			/// Divide a dual number by a real number
			inline dual2 operator/(real r) const {

				if(r == 0) {
					TH_MATH_ERROR("dual2::operator/", r, DIV_BY_ZERO);
					return dual2(nan(), nan(), nan());
				}

				return dual2(a / r, b / r, c / r);
			}


			/// Add a dual2 number from this one
			inline dual2& operator+=(const dual2& other) {

				a += other.a;
				b += other.b;
				c += other.c;
				return *this;
			}


			/// Sum a real number to this dual number
			inline dual2& operator+=(real r) {

				a += r;
				return *this;
			}

			/// Subtract a real number from this one
			inline dual2& operator-=(const dual2& other) {

				a -= other.a;
				b -= other.b;
				c -= other.c;
				return *this;
			}

			/// Subtract a real number from this dual number
			inline dual2& operator-=(real r) {

				a -= r;
				return *this;
			}

			/// Multiply this dual number by another one
			inline dual2& operator*=(const dual2& other) {

				a = (a * other.a);
				b = (a * other.b) + (b * other.a);
				c = (a * other.c) + (2 * b * other.b) + (c * other.a);
				return *this;
			}

			/// Multiply this dual number by a real number
			inline dual2& operator*=(real r) {
				a *= r;
				b *= r;
				c *= r;
				return *this;
			}

			/// Divide this dual number by another one
			inline dual2& operator/=(const dual2& other) {
				*this = operator*(other.inverse());
				return *this;
			}

			/// Divide a dual number by a real number
			inline dual2& operator/=(real r) {

				if(r == 0) {
					TH_MATH_ERROR("dual::operator/=", 0, DIV_BY_ZERO);
					a = nan();
					b = nan();
					c = nan();
					return *this;
				}

				a /= r;
				b /= r;
				c /= r;

				return *this;
			}


			/// Check whether two dual numbers have the same
			/// real and dual parts
			inline bool operator==(const dual2& other) {
				return (a == other.a) && (b == other.b) && (c == other.c);
			}


			/// Convert a dual number to a vector
			inline vec3 to_vec() const {
				vec3 res;
				res[0] = a;
				res[1] = b;
				res[2] = c;
				return res;
			}

			/// Initialize from a vector
			inline void from_vec(const vec3& v) {
				a = v[0];
				b = v[1];
				c = v[2];
			}


			// Friend operators to enable equations of the form
			// (real) op. (dual2)

			inline friend dual2 operator+(real a, const dual2& d) {
				return d + a;
			}

			inline friend dual2 operator-(real a, const dual2& d) {
				return -d + a;
			}

			inline friend dual2 operator*(real a, const dual2& d) {
				return d * a;
			}

			inline friend dual2 operator/(real a, const dual2& d) {
				return dual2(a, 0, 0) / d;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the dual number to string representation
			/// @param epsilon1 The character to use to represent epsilon1
			/// @param epsilon2 The character to use to represent epsilon2
			inline std::string to_string(const std::string& epsilon1 = "e1",
				const std::string& epsilon2 = "e2") const {

				std::stringstream res;

				res << a;
				res << (b >= 0 ? " + " : " - ");
				res << abs(b);

				res << epsilon1;

				res << (c >= 0 ? " + " : " - ");
				res << abs(c);

				res << epsilon2;

				return res.str();
			}


			/// Convert the dual number to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the dual number in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const dual2& obj) {
				return out << obj.to_string();
			}

#endif

	};

}


#endif
