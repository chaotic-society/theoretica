
///
/// @file dual.h Dual number class
///

#ifndef THEORETICA_DUAL_H
#define THEORETICA_DUAL_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/error.h"
#include "../core/constants.h"
#include "../algebra/algebra_types.h"


namespace theoretica {


	///
	/// @class dual
	/// Dual number class.
	/// Implemented as \f$a + b \epsilon\f$
	/// where \f$\epsilon\f$ is such that \f$\epsilon^2 = 0\f$
	///
	class dual {
		public:

			real a; // Real part
			real b; // "Dual" part

			/// Default constructor, initialize with null values
			dual() : a(0), b(0) {}

			/// Initialize from two real numbers
			dual(real real_part, real dual_part)
				: a(real_part), b(dual_part) {}

			/// Initialize from a real number
			dual(real real_part)
				: a(real_part), b(0) {}

			~dual() = default;
		
			/// Initialize from a vec2
			dual(const vec2& v) {
				a = v[0];
				b = v[1];
			}

			/// Initialize a dual number from a vec2
			inline dual& operator=(const vec2& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			/// Initialize a dual number from a real number
			inline dual& operator=(real x) {
				a = x;
				b = 0;
				return *this;
			}

			/// Initialize a dual number from a std::array
			inline dual& operator=(const std::array<real, 2>& v) {
				a = v[0];
				b = v[1];
				return *this;
			}

			/// Return real part
			inline const real& Re() const {
				return a;
			}

			/// Return real part
			inline real& Re() {
				return a;
			}

			/// Extract the real part of the dual number
			inline friend real Re(const dual& d) {
				return d.a;
			}


			/// Extract the real part of the dual number
			inline friend real& Re(dual& d) {
				return d.a;
			}

			/// Return dual part
			inline const real& Dual() const {
				return b;
			}

			/// Return dual part
			inline real& Dual() {
				return b;
			}

			/// Extract the real part of the dual number
			inline friend real Dual(const dual& d) {
				return d.b;
			}


			/// Extract the real part of the dual number
			inline friend real& Dual(dual& d) {
				return d.b;
			}

			/// Get the dual conjugate
			inline dual conjugate() const {
				return dual(a, -b);
			}

			/// Get the inverse of a dual number
			inline dual inverse() const {

				if(a == 0) {
					TH_MATH_ERROR("dual::inverse", 0, DIV_BY_ZERO);
					return dual(nan(), nan());
				}

				return dual(1.0 / a, -b / square(a));
			}

			/// Identity (for consistency)
			inline dual operator+() const {
				return dual(a, b);
			}

			/// Sum two dual numbers
			inline dual operator+(const dual& other) const {
				return dual(a + other.a, b + other.b);
			}

			/// Sum a real number to a dual number
			inline dual operator+(real r) const {
				return dual(a + r, b);
			}

			/// Get the opposite of a dual number
			inline dual operator-() const {
				return dual(-a, -b);
			}

			/// Subtract two dual numbers
			inline dual operator-(const dual& other) const {
				return dual(a - other.a, b - other.b);
			}

			/// Subtract a real number from a dual number
			inline dual operator-(real r) const {
				return dual(a - r, b);
			}

			/// Multiply two dual numbers
			inline dual operator*(const dual& other) const {
				return dual(a * other.a, a * other.b + b * other.a);
			}

			/// Multiply a dual number by a real number
			inline dual operator*(real r) const {
				return dual(a * r, b * r);
			}

			/// Dual division
			inline dual operator/(const dual& other) const {
				return dual(a / other.a,
					(b * other.a - a * other.b) / square(other.a));
			}

			/// Divide a dual number by a real number
			inline dual operator/(real r) const {
				return dual(a / r, b / r);
			}


			/// Add a dual number to this one 
			inline dual& operator+=(const dual& other) {

				a += other.a;
				b += other.b;
				return *this;
			}


			/// Sum a real number to this dual number
			inline dual& operator+=(real r) {

				a += r;
				return *this;
			}

			/// Subtract a real number from this one
			inline dual& operator-=(const dual& other) {

				a -= other.a;
				b -= other.b;
				return *this;
			}

			/// Subtract a real number from this dual number
			inline dual& operator-=(real r) {

				a -= r;
				return *this;
			}

			/// Multiply this dual number by another one
			inline dual& operator*=(const dual& other) {

				a = (a * other.a);
				b = (a * other.b) + (b * other.a);
				return *this;
			}

			/// Multiply this dual number by a real number
			inline dual& operator*=(real r) {

				a *= r;
				b *= r;
				return *this;
			}

			/// Divide this dual number by another one
			inline dual& operator/=(const dual& other) {

				a = (a / other.a);
				b = (b * other.a - a * other.b) / square(other.a);
				return *this;
			}

			/// Divide a dual number by a real number
			inline dual& operator/=(real r) {

				if(r == 0) {
					TH_MATH_ERROR("dual::operator/=", 0, DIV_BY_ZERO);
					a = nan();
					b = nan();
					return *this;
				}

				a /= r;
				b /= r;

				return *this;
			}


			/// Check whether two dual numbers have the same
			/// real and dual parts
			inline bool operator==(const dual& other) {
				return (a == other.a) && (b == other.b);
			}


			/// Convert a dual number to a vector
			inline vec2 to_vec() const {
				vec2 res;
				res[0] = a;
				res[1] = b;
				return res;
			}

			/// Initialize from a vector
			inline void from_vec(const vec2& v) {
				a = v[0];
				b = v[1];
			}


			/// Convert a dual number to matrix form
			inline mat2 to_mat() const {

				mat2 m;
				m(0, 0) = a;
				m(1, 0) = 0;
				m(0, 1) = b;
				m(1, 1) = a;
				return m;
			}


			// Friend operators to enable equations of the form
			// (real) op. (dual)

			inline friend dual operator+(real a, const dual& d) {
				return d + a;
			}

			inline friend dual operator-(real a, const dual& d) {
				return -d + a;
			}

			inline friend dual operator*(real a, const dual& d) {
				return d * a;
			}

			inline friend dual operator/(real a, const dual& d) {
				return dual(a, 0) / d;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the dual number to string representation
			/// @param epsilon The character to use to represent epsilon
			inline std::string to_string(const std::string& epsilon = "e") const {

				std::stringstream res;

				res << a;
				res << (b >= 0 ? " + " : " - ");
				res << abs(b) << epsilon;

				return res.str();
			}


			/// Convert the dual number to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the dual number in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const dual& obj) {
				return out << obj.to_string();
			}

#endif

	};

}


#endif
