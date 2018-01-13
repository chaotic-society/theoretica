#ifndef UROBORO_VEC_H
#define UROBORO_VEC_H
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

			inline vec4 operator*(real scalar) {
				return vec4(x * scalar, y * scalar, z * scalar, w * scalar);
			}

			inline void operator*=(real scalar) {
				x *= scalar;
				y *= scalar;
				z *= scalar;
			}

			inline void operator=(vec4 const& other) {
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


		class vec3 {

		public:
			real x, y, z;

			inline vec3() : x(0), y(0), z(0) {}

			inline vec3(real x, real y, real z)
			: x(x), y(y), z(z) {}

			inline ~vec3() {}

			inline real operator*(vec3 other) {
				return (x * other.x) + (y * other.y) + (z * other.z);
			}

			inline void operator*=(real scalar) {
				x *= scalar;
				y *= scalar;
				z *= scalar;
			}

			inline vec3 operator*(real scalar) {
				return vec3(x * scalar, y * scalar, z * scalar);
			}

			inline void operator=(vec3 const& other) {
				x = other.x;
				y = other.y;
				z = other.z;
			}

			inline bool operator==(vec3 other) {
				return (x == other.x) && (y == other.y)
					&& (z == other.z);
			}

			inline vec3 operator-(vec3 other) {
				return vec3(x - other.x, y - other.y, z - other.z);
			}

			inline void operator-=(vec3 other) {
				x -= other.x;
				y -= other.y;
				z -= other.z;
			}

			inline vec3 operator+(vec3 other) {
				return vec3(x + other.x, y + other.y, z + other.z);
			}

			inline void operator+=(vec3 other) {
				x += other.x;
				y += other.y;
				z += other.z;
			}

			inline bool operator<(vec3 other) {
				return sqrmagnitude() < other.sqrmagnitude();
			}

			inline bool operator>(vec3 other) {
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

			inline vec3 normalized() {
				real mag = magnitude();
				return vec3(x / mag, y / mag, z / mag);
			}

			inline vec3 cross(vec3 other) {
				return vec3(
					y * other.z - z * other.y,
					z * other.x - x * other.z,
					x * other.y - y * other.x
				);
			}

			inline real dot(vec3 other) {
				return (other.x * x) + (other.y * y) + (other.z * z);
			}

		};

		class vec2 {
			public:

				real x;
				real y;

				inline vec2(real x, real y) : x(x), y(y) {}

				inline vec2() : x(0), y(0) {}

				inline real operator*(vec2 other) {
					return (x * other.x) + (y * other.y);
				}

				inline vec2 operator*(real scalar) {
					return vec2(x * scalar, y * scalar);
				}

				inline void operator*=(real scalar) {
					x *= scalar;
					y *= scalar;
				}

				inline void operator=(vec2 const& other) {
					x = other.x;
					y = other.y;
				}

				inline bool operator==(vec2 other) {
					return (x == other.x) && (y == other.y);
				}

				inline vec2 operator-(vec2 other) {
					return vec2(x - other.x, y - other.y);
				}

				inline void operator-=(vec2 other) {
					x -= other.x;
					y -= other.y;
				}

				inline vec2 operator+(vec2 other) {
					return vec2(x + other.x, y + other.y);
				}

				inline void operator+=(vec2 other) {
					x += other.x;
					y += other.y;
				}

				inline bool operator<(vec2 other) {
					return sqrmagnitude() < other.sqrmagnitude();
				}

				inline bool operator>(vec2 other) {
					return sqrmagnitude() > other.sqrmagnitude();
				}

				inline real magnitude() {
					return sqrt((x * x) + (y * y));
				}

				inline real sqrmagnitude() {
					return (x * x) + (y * y);
				}

				inline void normalize() {
					real mag = magnitude();
					x /= mag;
					y /= mag;
				}

				inline vec2 normalized() {
					real mag = magnitude();
					return vec2(x / mag, y / mag);
				}

				inline real cross(vec2 other) {
					return x * other.y - y * other.x;
				}

				inline real dot(vec2 other) {
					return (other.x * x) + (other.y * y);
				}

		};


}

#endif
