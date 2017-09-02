#ifndef _VEC_H
#define _VEC_H
#include "common.h"

namespace uroboro {

	class vec4 {

		public:
			real x, y, z, w;

			inline vec4() : x(0), y(0), z(0), w(0) {}

			inline vec4(real x, real y, real z, real w = 0)
			: x(x), y(y), z(z), w(w) {}

			inline ~vec4() {}

			inline real operator*(vec4 other) {
				return (x * other.x) + (y * other.y) + (z * other.z);
			}

			inline void operator*=(real scalar) {
				x *= scalar;
				y *= scalar;
				z *= scalar;
			}

			inline void operator=(vec4 &other) {
				x = other.x;
				y = other.y;
				z = other.z;
				w = other.w;
			}

			inline bool operator==(vec4 other) {
				return (x == other.x) && (y == other.y)
					&& (z == other.z) && (w == other.w);
			}

			inline vec4 operator-(vec4 other) {
				return vec4(x - other.x, y - other.y, z - other.z, w);
			}

			inline void operator-=(vec4 other) {
				x -= other.x;
				y -= other.y;
				z -= other.z;
			}

			inline vec4 operator+(vec4 other) {
				return vec4(x + other.x, y + other.y, z + other.z, w);
			}

			inline void operator+=(vec4 other) {
				x += other.x;
				y += other.y;
				z += other.z;
			}

			inline bool operator<(vec4 other) {
				return sqrmagnitude() < other.sqrmagnitude();
			}

			inline bool operator>(vec4 other) {
				return sqrmagnitude() > other.sqrmagnitude();
			}

			inline real magnitude() {
				return sqrt((x * x) + (y * y) + (z * z));
			}

			inline real sqrmagnitude() {
				return (x * x) + (y * y) + (z * z);
			}

			inline void normalize() {
				real mag = magnitude();
				x /= mag;
				y /= mag;
				z /= mag;
			}

			inline vec4 normalized() {
				real mag = magnitude();
				return vec4(x / mag, y / mag, z / mag);
			}

			inline vec4 cross(vec4 other) {
				return vec4(
					y * other.z - z * other.y,
					z * other.x - x * other.z,
					x * other.y - y * other.x,
					w
				);
			}

			inline real dot(vec4 other) {
				return (other.x * x) + (other.y * y) + (other.z * z);
			}

		};


}

#endif
