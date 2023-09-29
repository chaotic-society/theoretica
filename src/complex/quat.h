
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
	/// Quaternion class in the form \f$a + bi + cj + dk\f$
	template<typename Type = real>
	class quat {
		public:

			/// Real part
			Type a;

			/// Imaginary parts
			Type b, c, d;


			/// Construct a quaternion with all zero components
			quat() : a(0), b(0), c(0), d(0) {}
			

			/// Construct a quaternion from a real number
			quat(Type r) : a(r), b(0), c(0), d(0) {}


			/// Construct a quaternion from its components
			quat(Type a, Type b, Type c, Type d)
			: a(a), b(b), c(c), d(d) {}


			/// Copy constructor
			quat(const quat& q) : a(q.a), b(q.b), c(q.c), d(q.d) {}


			/// Assignment operator
			inline quat& operator=(const quat& q) {
				a = q.a;
				b = q.b;
				c = q.c;
				d = q.d;
				return *this;
			}


			/// Assignment operator from a 4D array
			template<typename T>
			inline quat& operator=(const std::array<T, 4>& v) {
				a = v[0];
				b = v[1];
				c = v[2];
				d = v[3];
				return *this;
			}


			/// Get the real part of the quaternion
			inline Type Re() const {
				return a;
			}


			/// Extract the real part of the quaternion
			inline friend Type Re(const quat& q) {
				return q.a;
			}


			/// Get the first imaginary part of the quaternion
			inline Type Im1() const {
				return b;
			}


			/// Extract the first imaginary part of the quaternion
			inline friend Type Im1(const quat& q) {
				return q.b;
			}


			/// Get the second imaginary part of the quaternion
			inline Type Im2() const {
				return c;
			}


			/// Extract the second imaginary part of the quaternion
			inline friend Type Im2(const quat& q) {
				return q.c;
			}


			/// Get the third imaginary part of the quaternion
			inline Type Im3() const {
				return d;
			}


			/// Extract the third imaginary part of the quaternion
			inline friend Type Im3(const quat& q) {
				return q.d;
			}


			/// Compute the conjugate of the quaternion
			inline quat conjugate() const {
				return quat(a, -b, -c, -d);
			}


			/// Compute the square norm of the quaternion
			inline Type sqr_norm() const {
				return a * a + b * b + c * c + d * d;
			}


			/// Compute the norm of the quaternion
			inline Type norm() const {
				return sqrt(sqr_norm());
			}


			/// Compute the inverse of the quaternion
			inline quat inverse() const {

				const Type n = sqr_norm();

				if(n < MACH_EPSILON) {
					TH_MATH_ERROR("quat::inverse", n, DIV_BY_ZERO);
					return quat((Type) nan());
				}

				return conjugate() / sqr_norm;
			}


			/// Invert the quaternion
			inline quat& invert() {
				
				const Type n = sqr_norm();

				if(n < MACH_EPSILON) {
					TH_MATH_ERROR("quat::invert", n, DIV_BY_ZERO);
					return quat((Type) nan());
				}

				return (*this = quat(
					a / n,
					b / n * -1,
					c / n * -1,
					d / n * -1
				));
			}


			/// Identity (for consistency)
			inline quat operator+() const {
				return *this;
			}


			/// Get the opposite of the quaternion
			inline quat operator-() const {
				return quat(-a, -b, -c, -d);
			}


			/// Sum a quaternion and a scalar value
			inline quat operator+(Type k) const {
				return quat(a + k, b, c, d);
			}


			/// Subtract a quaternion and a scalar value
			inline quat operator-(Type k) const {
				return quat(a - k, b, c, d);
			}


			/// Multiply a quaternion by a scalar value
			inline quat operator*(Type k) const {
				return quat(a * k, b * k, c * k, d * k);
			}


			/// Divide a quaternion by a scalar value
			inline quat operator/(Type k) const {

				if(abs(k) < MACH_EPSILON) {
					TH_MATH_ERROR("quat::operator/", k, DIV_BY_ZERO);
					return quat(nan());
				}

				return quat(a / k, b / k, c / k, d / k);
			}


			/// Add two quaternions
			inline quat operator+(const quat& other) const {
				return quat(a + other.a, b + other.b, c + other.c, d + other.d);
			}


			/// Subtract two quaternions
			inline quat operator-(const quat& other) const {
				return quat(a - other.a, b - other.b, c - other.c, d - other.d);
			}


			/// Multiply two quaternions
			inline quat operator*(const quat& q) const {

				quat r;

				r.a = a * q.a - b * q.b - c * q.c - d * q.d;
				r.b = a * q.b + b * q.a + c * q.d - d * q.c;
				r.c = a * q.c - b * q.d + c * q.a + d * q.b;
				r.d = a * q.d + b * q.c - c * q.b + d * q.a;

				return r;
			}


			/// Divide two quaternions
			inline quat operator/(const quat& other) const {
				return operator*(other.inverse());
			}


			/// Multiply this quaternion by a scalar value
			inline quat& operator*=(Type k) {
				return (*this = operator*(k));
			}


			/// Divide this quaternion by a scalar value
			inline quat& operator/=(Type k) {

				if(abs(k) < MACH_EPSILON) {
					TH_MATH_ERROR("quat::operator/=", k, DIV_BY_ZERO);
					return (*this = quat(nan()));
				}

				a /= k;
				b /= k;
				c /= k;
				d /= k;
				return *this;
			}


			/// Add a quaternion to this one
			inline quat& operator+=(const quat& other) {
				a += other.a;
				b += other.b;
				c += other.c;
				d += other.d;
				return *this;
			}


			/// Subtract a quaternion to this one
			inline quat& operator-=(const quat& other) {
				a -= other.a;
				b -= other.b;
				c -= other.c;
				d -= other.d;
				return *this;
			}


			/// Multiply this quaternion by another one
			inline quat& operator*=(const quat& other) {
				return (*this = operator*(other));
			}


			/// Divide this quaternion by another one
			inline quat& operator/=(const quat& other) {
				return operator*=(other.inverse());
			}


			/// Transform a 3D vector
			template<typename Vector>
			inline Vector transform(const Vector& v) const {

				Vector res;
				res.resize(3);

				if(v.size() != 3) {
					TH_MATH_ERROR("quat::transform", v.size(), INVALID_ARGUMENT);
					vec_error(res);
					return res;
				}

				const quat q = quat(0, v.get(0), v.get(1), v.get(2));
				const quat r = (*this * q) * inverse();

				res[0] = r.b;
				res[1] = r.c;
				res[2] = r.d;

				return res;
			}


			/// Construct a quaternion which represents a rotation
			/// of <rad> radians around the <axis> arbitrary axis
			template<typename Vector>
			inline static quat rotation(Type rad, const Vector& axis) {
				return quat(
					cos(rad / 2.0),
					algebra::normalize(axis) * sin(rad / 2.0)
				);
			}


			/// Rotate a 3D vector <v> by <rad> radians around
			/// the <axis> arbitrary axis
			template<typename Vector>
			inline static Vector rotate(const Vector& v, Type rad, const Vector& axis) {

				Vector res;
				res.resize(3);

				if(axis.size() != 3) {
					TH_MATH_ERROR("quat::rotate", axis.size(), INVALID_ARGUMENT);
					vec_error(res);
					return res;
				}

				if(v.size() != 3) {
					TH_MATH_ERROR("quat::rotate", v.size(), INVALID_ARGUMENT);
					vec_error(res);
					return res;
				}

				Vector n_axis = algebra::normalize(axis);

				const Type s = sin(rad / 2.0);
				const Type c = cos(rad / 2.0);

				const quat q = quat(c, n_axis[0] * s, n_axis[1] * s, n_axis[2] * s);
				const quat q_inv = q.conjugate();
				const quat p = quat(0, v[0], v[1], v[2]);

				const quat r = q * p * q_inv;
				res[0] = r.b;
				res[1] = r.c;
				res[2] = r.d;

				return res;
			}


			// Friend operators to enable equations of the form
			// (real) op. (quat)

			inline friend quat operator+(Type r, const quat& z) {
				return z + quat(r, 0, 0, 0);
			}

			inline friend quat operator-(Type r, const quat& z) {
				return (z * -1) + quat(r, 0, 0, 0);
			}

			inline friend quat operator*(Type r, const quat& z) {
				return z * r;
			}

			inline friend quat operator/(Type r, const quat& z) {
				return quat(r, 0, 0, 0) / z;
			}




#ifndef THEORETICA_NO_PRINT

			/// Convert the quaternion to string representation
			inline std::string to_string() const {

				std::stringstream res;

				res << a;
				res << (b >= 0 ? " + " : " - ") << abs(b) << "i";
				res << (c >= 0 ? " + " : " - ") << abs(c) << "j";
				res << (d >= 0 ? " + " : " - ") << abs(d) << "k";

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
