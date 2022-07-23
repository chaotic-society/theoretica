
///
/// @file quat.h Quaternion algebra
///

#ifndef THEORETICA_QUATERNION_H
#define THEORETICA_QUATERNION_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/real_analysis.h"
#include "../algebra/vec.h"


namespace theoretica {

	/// @class quat
	/// Quaternion in the form \f$a + bi + cj + dk\f$,
	/// stored as \f$(a + \vec v)\f$ where \f$a \in \mathbb{R}\f$ and
	/// \f$\vec v \in \mathbb{R}^3\f$
	class quat {
		public:

			real a;
			vec3 v;

			/// Initialize as (0 + 0i + 0j + 0k)
			quat() : a(0), v(vec3()) {}

			/// Initialize from a real number and a vector
			quat(real a, const vec3& v) : a(a), v(v) {}

			/// Initialize from another quaternion
			quat(const quat& other) : a(other.a), v(other.v) {}

			inline quat& operator=(const quat& other) {
				a = other.a;
				v = other.v;
				return *this;
			}

			inline quat& operator=(const std::array<real, 4>& v) {
				a = v[0];
				this->v.data[0] = v[0];
				this->v.data[1] = v[1];
				this->v.data[2] = v[2];
				return *this;
			}

			/// Initialize from four real numbers
			quat(real a, real b, real c, real d) : a(a) {
				v.data[0] = b;
				v.data[1] = c;
				v.data[2] = d;
			}

			~quat() = default;

			/// Get the norm of a quaternion
			inline real norm() const {
				return sqrt(a * a + v.square_magnitude());
			}

			/// Get the square norm of a quaternion
			inline real square_norm() const {
				return a * a + v.square_magnitude();
			}

			/// Return the conjugate of a quaternion
			inline quat conjugate() const {
				return quat(-a, v * -1);
			}

			// Operators

			/// Multiply a quaternion by a scalar value
			inline quat operator*(real scalar) const {
				return quat(a * scalar, v * scalar);
			}

			/// Divide a quaternion by a scalar value
			inline quat operator/(real scalar) const {

				if(scalar == 0) {
					TH_MATH_ERROR("quat::operator/", scalar, DIV_BY_ZERO);
					return quat(nan(), vec3(nan()));
				}

				return quat(a / scalar, v / scalar);
			}

			/// Add two quaternions
			inline quat operator+(const quat& other) const {
				return quat(a + other.a, v + other.v);
			}

			/// Subtract two quaternions
			inline quat operator-(const quat& other) const {
				return quat(a - other.a, v - other.v);
			}

			/// Multiply two quaternions
			inline quat operator*(const quat& other) const {
				return quat((a * other.a) - (v * other.v),
							(other.v * a) + (v * other.a) + v.cross(other.v));
			}

			/// Divide two quaternions
			inline quat operator/(const quat& other) const {
				return operator*(other.inverse());
			}


			/// Multiply this quaternion by a scalar value
			inline quat& operator*=(real scalar) {

				a *= scalar;
				v *= scalar;
				return *this;
			}

			/// Divide this quaternion by a scalar value
			inline quat& operator/=(real scalar) {

				if(scalar == 0) {
					TH_MATH_ERROR("quat::operator/=", scalar, DIV_BY_ZERO);
					return (*this = quat(nan(), vec3(nan())));
				}

				a /= scalar;
				v /= scalar;
				return *this;
			}

			/// Add a quaternion to this one
			inline quat& operator+=(const quat& other) {

				a += other.a;
				v += other.v;
				return *this;
			}

			/// Subtract a quaternion to this one
			inline quat& operator-=(const quat& other) {

				a -= other.a;
				v -= other.v;
				return *this;
			}

			/// Multiply this quaternion by another one
			inline quat& operator*=(const quat& other) {

				a = (a * other.a) - (v * other.v);
				v = (other.v * a) + (v * other.a) + v.cross(other.v);
				return *this;
			}

			/// Divide this quaternion by another one
			inline quat& operator/=(const quat& other) {
				return operator*=(other.inverse());
			}


			/// Normalize the quaternion
			inline void normalize() {

				real n = norm();

				if(n == 0) {
					TH_MATH_ERROR("quat::normalize", n, DIV_BY_ZERO);
					*this = quat(nan(), vec3(nan()));
				}

				a /= n;
				v /= n;
			}

			/// Return the normalized quaternion
			inline quat normalized() const {

				real n = norm();

				if(n == 0) {
					TH_MATH_ERROR("quat::normalized", n, DIV_BY_ZERO);
					return quat(nan(), vec3(nan()));
				}

				return quat(a / n, v / n);
			}

			/// Return the inverse of a quaternion
			inline quat inverse() const {

				real sqr_norm = square_norm();

				if(sqr_norm == 0) {
					TH_MATH_ERROR("quat::inverse", sqr_norm, DIV_BY_ZERO);
					return quat(nan(), vec3(nan()));
				}

				return conjugate() / sqr_norm;
			}


			/// Convert a quaternion to a vector in the form
			/// (a, b, c, d)
			inline vec4 to_vec4() const {
				return vec4({a, v.get(0), v.get(1), v.get(2)});
			}

			/// Convert the quaternion to a 4x4 matrix
			inline mat4 to_mat4() const {

				real x = v.get(0);
				real y = v.get(1);
				real z = v.get(2);
				real w = a;

				mat4 res;

				res.iat(0, 0) = 1 - (2 * square(y) + 2 * square(z));
				res.iat(0, 1) = 2 * x * y + 2 * z * w;
				res.iat(0, 2) = 2 * x * z - 2 * y * w;
				res.iat(0, 3) = 0;

				res.iat(1, 0) = 2 * x * y - 2 * z * w;
				res.iat(1, 1) = 1 - (2 * square(x) + 2 * square(z));
				res.iat(1, 2) = 2 * y * z + 2 * x * w;
				res.iat(1, 3) = 0;

				res.iat(2, 0) = 2 * x * z + 2 * y * w;
				res.iat(2, 1) = 2 * y * z - 2 * x * w;
				res.iat(2, 2) = 1 - (2 * square(x) + 2 * square(y));
				res.iat(2, 3) = 0;

				res.iat(3, 0) = 0;
				res.iat(3, 1) = 0;
				res.iat(3, 2) = 0;
				res.iat(3, 3) = 1;

				return res;
			}

			/// Convert the quaternion to a 3x3 matrix
			inline mat3 to_mat3() const {

				real x = v.get(0);
				real y = v.get(1);
				real z = v.get(2);
				real w = a;

				mat3 res;

				res.iat(0, 0) = 1 - (2 * square(y) + 2 * square(z));
				res.iat(0, 1) = 2 * x * y + 2 * z * w;
				res.iat(0, 2) = 2 * x * z - 2 * y * w;

				res.iat(1, 0) = 2 * x * y - 2 * z * w;
				res.iat(1, 1) = 1 - (2 * square(x) + 2 * square(z));
				res.iat(1, 2) = 2 * y * z + 2 * x * w;

				res.iat(2, 0) = 2 * x * z + 2 * y * w;
				res.iat(2, 1) = 2 * y * z - 2 * x * w;
				res.iat(2, 2) = 1 - (2 * square(x) + 2 * square(y));

				return res;
			}


			/// Transform a 3D vector
			inline vec3 transform(const vec3& v) const {
				return (operator*(quat(0, v)) * inverse()).v;
			}


			/// Construct a quaternion which represents a rotation
			/// of <rad> radians around the <axis> arbitrary axis
			inline static quat rotation(real rad, const vec3& axis) {
				return quat(cos(rad / 2.0),
							axis.normalized() * sin(rad / 2.0));
			}


			/// Rotate a 3D vector <v> by <rad> radians around
			/// the <axis> arbitrary axis
			inline static vec3 rotate(const vec3& v, real rad, const vec3& axis) {

				vec3 n_axis = axis.normalized();

				quat q = quat(cos(rad / 2.0),
					n_axis * sin(rad / 2.0));

				quat q_inv = quat(cos(rad / 2.0),
					n_axis * (sin(rad / 2.0) * -1));

				quat p = quat(0, v);

				return (q * p * q_inv).v;
			}


			// Friend operators to enable equations of the form
			// (real) op. (quat)

			inline friend quat operator+(real r, const quat& z) {
				return z + quat(r, 0, 0, 0);
			}

			inline friend quat operator-(real r, const quat& z) {
				return (z * -1) + quat(r, 0, 0, 0);
			}

			inline friend quat operator*(real r, const quat& z) {
				return z * r;
			}

			inline friend quat operator/(real r, const quat& z) {
				return quat(r, 0, 0, 0) / z;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the quaternion to string representation
			inline std::string to_string() const {

				std::stringstream res;

				res << a;
				res << (v.get(0) >= 0 ? " + " : " - ") << "i" << abs(v.get(0));
				res << (v.get(1) >= 0 ? " + " : " - ") << "j" << abs(v.get(1));
				res << (v.get(2) >= 0 ? " + " : " - ") << "k" << abs(v.get(2));

				return res.str();
			}


			/// Stream the quaternion in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const quat& obj) {
				return out << obj.to_string();
			}

#endif

	};

}

#endif
