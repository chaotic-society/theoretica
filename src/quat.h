#ifndef UROBORO_QUATERNION_H
#define UROBORO_QUATERNION_H

#include "./common.h"
#include "./vec.h"

namespace uroboro {

	// Quaternion implementation
	// In the form a + bi + cj + dk
	// As a + v [real + vec<3>]
	class quat {
		public:

			real a;
			vec3 v;

			inline quat() {
				a = 0;
				v = vec3();
			}

			inline quat(real a, const vec3& v) {
				this->a = a;
				this->v = v;
			}

			inline quat(const quat& other) {
				a = other.a;
				v = other.v;
			}

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

			inline quat(real a, real b, real c, real d) {
				this-> a = a;
				v.data[0] = b;
				v.data[1] = c;
				v.data[2] = d;
			}

			inline ~quat() {}

			inline real norm() const {
				return uroboro::sqrt(a * a + v.square_magnitude());
			}

			inline real square_norm() const {
				return a * a + v.square_magnitude();
			}

			inline quat conjugate() const {
				return quat(-a, v * -1);
			}

			// Operators

			inline quat operator*(real scalar) {
				return quat(a * scalar, v * scalar);
			}

			inline quat operator/(real scalar) {
				return quat(a / scalar, v / scalar);
			}

			inline quat operator+(const quat& other) const {
				return quat(a + other.a, v + other.v);
			}

			inline quat operator-(const quat& other) const {
				return quat(a - other.a, v - other.v);
			}

			inline quat operator*(const quat& other) const {
				return quat((a * other.a) - (v * other.v),
							(other.v * a) + (v * other.a) + v.cross(other.v));
			}

			inline quat operator/(const quat& other) {
				return operator*(other.inverse());
			}

			inline void normalize() {
				real n = norm();
				a /= n;
				v = v / n;
			}

			inline quat normalized() const {
				real n = norm();
				return quat(a / n, v / n);
			}

			inline quat inverse() const {
				return conjugate() / square_norm();
			}

			// Transform a 3D vector
			inline vec3 transform(const vec3& v) const {
				return (operator*(quat(0, v)) * inverse()).v;
			}

			// Return a quaternion which represents a rotation
			// of <rad> radians around the <axis> arbitrary axis
			inline static quat rotation(real rad, const vec3& axis) {
				return quat(uroboro::cos(rad / 2.0),
							axis.normalized() * uroboro::sin(rad / 2.0));
			}

			// Rotate a 3D vector <v> by <rad> radians around
			// the <axis> arbitrary axis
			inline static vec3 rotate(const vec3& v, real rad, const vec3& axis) {

				vec3 n_axis = axis.normalized();

				quat q = quat(uroboro::cos(rad / 2.0),
					n_axis * uroboro::sin(rad / 2.0));

				quat q_inv = quat(uroboro::cos(rad / 2.0),
					n_axis * (uroboro::sin(rad / 2.0) * -1));

				quat p = quat(0, v);

				return (q * p * q_inv).v;
			}

	};

}

#endif
