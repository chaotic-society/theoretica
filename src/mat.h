#ifndef UROBORO_MAT_H
#define UROBORO_MAT_H
#include "vec.h"

namespace uroboro {

	class mat4 {

		public:
			real data[4][4] = { {0.f}, {0.f}, {0.f}, {0.f} }; //4x4

			inline mat4() {
				identity();
			}

			inline mat4(real i) {
				data[0][0] = i;
				data[1][1] = i;
				data[2][2] = i;
				data[3][3] = i;
			}

			inline mat4(real a, real b, real c, real d,
						real e, real f, real g, real h,
						real i, real j, real k, real l,
						real m, real n, real o, real p) {
				data[0][0] = a;
				data[1][0] = b;
				data[2][0] = c;
				data[3][0] = d;
				data[0][1] = e;
				data[1][1] = f;
				data[2][1] = g;
				data[3][1] = h;
				data[0][2] = i;
				data[1][2] = j;
				data[2][2] = k;
				data[3][2] = l;
				data[0][3] = m;
				data[1][3] = n;
				data[2][3] = o;
				data[3][3] = p;
			}

			inline ~mat4() {}

			inline void identity() {
				data[0][0] = 1;
				data[0][1] = 0;
				data[0][2] = 0;
				data[0][3] = 0;
				data[1][0] = 0;
				data[1][1] = 1;
				data[1][2] = 0;
				data[1][3] = 0;
				data[2][0] = 0;
				data[2][1] = 0;
				data[2][2] = 1;
				data[2][3] = 0;
				data[3][0] = 0;
				data[3][1] = 0;
				data[3][2] = 0;
				data[3][3] = 1;
			}

			inline void invert() {

			    const real det = data[0][2] * data[1][1] * data[2][0]
			                    + data[0][1] * data[1][2] * data[2][0]
			                    - data[0][2] * data[1][0] * data[2][1]
			                    - data[0][0] * data[1][2] * data[2][1]
			                    - data[0][1] * data[1][0] * data[2][2]
			                    + data[0][0] * data[1][1] * data[2][2];

			    if(det == 0.f) return;
			    const real invdet = 1.f/det;

			    data[0][0] = (-data[1][2] * data[2][1] + data[1][1] * data[2][2]) * invdet;
			    data[0][1] = (data[0][2] * data[2][1] - data[0][1] * data[2][2]) * invdet;
			    data[0][2] = (-data[0][2] * data[1][1] + data[0][1] * data[1][2] * data[3][3]) * invdet;
			    data[1][0] = (data[1][2] * data[2][0] - data[1][0] * data[2][2]) * invdet;
			    data[1][1] = (-data[0][2] * data[2][0] + data[0][0] * data[2][2]) * invdet;
			    data[1][2] = (data[0][2] * data[1][0] - data[0][0] * data[1][2] * data[3][3]) * invdet;
			    data[2][0] = (-data[1][1] * data[2][0] + data[1][0] * data[2][1] * data[3][3]) * invdet;
			    data[2][1] = (data[0][1] * data[2][0] - data[0][0] * data[2][1] * data[3][3]) * invdet;
			    data[2][2] = (-data[0][1] * data[1][0] + data[0][0] * data[1][1] * data[3][3]) * invdet;
			    data[3][0] = (data[1][2] * data[2][1] * data[3][0]
	                        - data[1][1] * data[2][2] * data[3][0]
	                        - data[1][2] * data[2][0] * data[3][1]
	                        + data[1][0] * data[2][2] * data[3][1]
	                        + data[1][1] * data[2][0] * data[3][2]
	                        - data[1][0] * data[2][1] * data[3][2]) * invdet;
			    data[3][1] = (data[0][2] * data[2][1] * data[3][0]
	                        - data[0][1] * data[2][2] * data[3][0]
	                        - data[0][2] * data[2][0] * data[3][1]
	                        + data[0][0] * data[2][2] * data[3][1]
	                        + data[0][1] * data[2][0] * data[3][2]
	                        - data[0][0] * data[2][1] * data[3][2]) * invdet;
			    data[3][2] = (data[0][2] * data[1][1] * data[3][0]
	                        - data[0][1] * data[1][2] * data[3][0]
	                        - data[0][2] * data[1][0] * data[3][1]
	                        + data[0][0] * data[1][2] * data[3][1]
	                        + data[0][1] * data[1][0] * data[3][2]
	                        - data[0][0] * data[1][1] * data[3][2]) * invdet;
			}

			inline void transpose() {
				real buffer[16];
				buffer[0] = data[0][0];
				buffer[1] = data[1][0];
				buffer[2] = data[2][0];
				buffer[3] = data[3][0];
				buffer[4] = data[0][1];
				buffer[5] = data[1][1];
				buffer[6] = data[2][1];
				buffer[7] = data[3][1];
				buffer[8] = data[0][2];
				buffer[9] = data[1][2];
				buffer[10] = data[2][2];
				buffer[11] = data[3][2];
				buffer[12] = data[0][3];
				buffer[13] = data[1][3];
				buffer[14] = data[2][3];
				buffer[15] = data[3][3];

				data[0][1] = buffer[1];
				data[0][2] = buffer[2];
				data[0][3] = buffer[3];
				data[1][0] = buffer[4];
				data[1][1] = buffer[5];
				data[1][2] = buffer[6];
				data[1][3] = buffer[7];
				data[2][0] = buffer[8];
				data[2][1] = buffer[9];
				data[2][3] = buffer[11];
				data[3][0] = buffer[12];
				data[3][1] = buffer[13];
				data[3][2] = buffer[14];
			}

			inline vec4 transform(vec4 vector) {
				return vec4((data[0][0] * vector.x) + (data[1][0] * vector.y) + (data[2][0] * vector.z) + (data[3][0] * vector.w),
							(data[0][1] * vector.x) + (data[1][1] * vector.y) + (data[2][1] * vector.z) + (data[3][1] * vector.w),
							(data[0][2] * vector.x) + (data[1][2] * vector.y) + (data[2][2] * vector.z) + (data[3][2] * vector.w));
			}

			inline vec3 transform(vec3 vector) {
				return vec3((data[0][0] * vector.x) + (data[1][0] * vector.y) + (data[2][0] * vector.z) + data[3][0],
							(data[0][1] * vector.x) + (data[1][1] * vector.y) + (data[2][1] * vector.z) + data[3][1],
							(data[0][2] * vector.x) + (data[1][2] * vector.y) + (data[2][2] * vector.z) + data[3][2]);
			}

			inline void translate(float x, float y, float z) {
				data[3][0] += x;
				data[3][1] += y;
				data[3][2] += z;
			}

			inline void translate(vec4 translation) {
				data[3][0] += translation.x;
				data[3][1] += translation.y;
				data[3][2] += translation.z;
			}

			inline void translate(vec3 translation) {
				data[3][0] += translation.x;
				data[3][1] += translation.y;
				data[3][2] += translation.z;
			}

			inline void rotate(real radians, vec3 rotation) {
			}

			//inline void rotate(quat rotation) {}

			inline void scale(float x, float y, float z) {
				data[0][0] *= x;
				data[1][1] *= y;
				data[2][2] *= z;
			}

			inline void scale(vec4 scale) {
				data[0][0] *= scale.x;
				data[1][1] *= scale.y;
				data[2][2] *= scale.z;
			}

			inline void scale(vec3 scale) {
				data[0][0] *= scale.x;
				data[1][1] *= scale.y;
				data[2][2] *= scale.z;
			}

			inline void operator=(const mat4 &other) {
				data[0][0] = other.data[0][0];
				data[0][1] = other.data[0][1];
				data[0][2] = other.data[0][2];
				data[0][3] = other.data[0][3];
				data[1][0] = other.data[1][0];
				data[1][1] = other.data[1][1];
				data[1][2] = other.data[1][2];
				data[1][3] = other.data[1][3];
				data[2][0] = other.data[2][0];
				data[2][1] = other.data[2][1];
				data[2][2] = other.data[2][2];
				data[2][3] = other.data[2][3];
				data[3][0] = other.data[3][0];
				data[3][1] = other.data[3][1];
				data[3][2] = other.data[3][2];
				data[3][3] = other.data[3][3];
			}

			//inline void operator=(const quat &quaternion) {}

			inline mat4 operator*(real scalar) {

				mat4 result;

				result.data[0][0] = data[0][0] * scalar;
				result.data[0][1] = data[0][1] * scalar;
				result.data[0][2] = data[0][2] * scalar;
				result.data[0][3] = data[0][3] * scalar;
				result.data[1][0] = data[1][0] * scalar;
				result.data[1][1] = data[1][1] * scalar;
				result.data[1][2] = data[1][2] * scalar;
				result.data[1][3] = data[1][3] * scalar;
				result.data[2][0] = data[2][0] * scalar;
				result.data[2][1] = data[2][1] * scalar;
				result.data[2][2] = data[2][2] * scalar;
				result.data[2][3] = data[2][3] * scalar;
				result.data[3][0] = data[3][0] * scalar;
				result.data[3][1] = data[3][1] * scalar;
				result.data[3][2] = data[3][2] * scalar;
				result.data[3][3] = data[3][3] * scalar;

				return result;
			}

			inline vec4 operator*(vec4 vector) {
				return transform(vector);
			}

			inline vec3 operator*(vec3 vector) {
				return transform(vector);
			}

			inline mat4 operator*(const mat4 &other) {
				return mat4(other.data[0][0] * data[0][0] + other.data[0][1] * data[1][0] + other.data[0][2] * data[2][0],
							other.data[1][0] * data[0][0] + other.data[1][1] * data[1][0] + other.data[1][2] * data[2][0],
							other.data[2][0] * data[0][0] + other.data[2][1] * data[1][0] + other.data[2][2] * data[2][0],
							other.data[3][0] * data[0][0] + other.data[3][1] * data[1][0] + other.data[3][2] * data[2][0] + data[3][0],
							other.data[0][0] * data[0][1] + other.data[0][1] * data[1][1] + other.data[0][2] * data[2][1],
							other.data[1][0] * data[0][1] + other.data[1][1] * data[1][1] + other.data[1][2] * data[2][1],
							other.data[2][0] * data[0][1] + other.data[2][1] * data[1][1] + other.data[2][2] * data[2][1],
							other.data[3][0] * data[0][1] + other.data[3][1] * data[1][1] + other.data[3][2] * data[2][1] + data[3][1],
							other.data[0][0] * data[0][2] + other.data[0][1] * data[1][2] + other.data[0][2] * data[2][2],
							other.data[1][0] * data[0][2] + other.data[1][1] * data[1][2] + other.data[1][2] * data[2][2],
							other.data[2][0] * data[0][2] + other.data[2][1] * data[1][2] + other.data[2][2] * data[2][2],
							other.data[3][0] * data[0][2] + other.data[3][1] * data[1][2] + other.data[3][2] * data[2][2] + data[3][2],
							0, 0, 0, 1
				);
			}

			inline mat4 operator+(real scalar) {

				mat4 result;

				result.data[0][0] = data[0][0] + scalar;
				result.data[0][1] = data[0][1] + scalar;
				result.data[0][2] = data[0][2] + scalar;
				result.data[0][3] = data[0][3] + scalar;
				result.data[1][0] = data[1][0] + scalar;
				result.data[1][1] = data[1][1] + scalar;
				result.data[1][2] = data[1][2] + scalar;
				result.data[1][3] = data[1][3] + scalar;
				result.data[2][0] = data[2][0] + scalar;
				result.data[2][1] = data[2][1] + scalar;
				result.data[2][2] = data[2][2] + scalar;
				result.data[2][3] = data[2][3] + scalar;
				result.data[3][0] = data[3][0] + scalar;
				result.data[3][1] = data[3][1] + scalar;
				result.data[3][2] = data[3][2] + scalar;
				result.data[3][3] = data[3][3] + scalar;

				return result;
			}

			// inline real& operator[](unsigned int i) {
			// 	return data[i];
			// }

	};


	inline mat4 translation(real x, real y, real z) {

		mat4 result = mat4();

		result.data[0][3] = x;
		result.data[1][3] = y;
		result.data[2][3] = z;

		return result;
	}

	inline mat4 translation(vec4 translation) {

		mat4 result = mat4();

		result.data[0][3] = translation.x;
		result.data[1][3] = translation.y;
		result.data[2][3] = translation.z;

		return result;
	}

	inline mat4 translation(vec3 translation) {

		mat4 result = mat4();

		result.data[0][3] = translation.x;
		result.data[1][3] = translation.y;
		result.data[2][3] = translation.z;

		return result;
	}

	//inline mat4 rotation() {}

	inline mat4 scale(real x, real y, real z) {

		mat4 result = mat4();

		result.data[0][0] = x;
		result.data[1][1] = y;
		result.data[2][2] = z;

		return result;
	}

	inline mat4 scale(vec4 scale) {

		mat4 result = mat4();

		result.data[0][0] = scale.x;
		result.data[1][1] = scale.y;
		result.data[2][2] = scale.z;

		return result;
	}

	inline mat4 scale(vec3 scale) {

		mat4 result = mat4();

		result.data[0][0] = scale.x;
		result.data[1][1] = scale.y;
		result.data[2][2] = scale.z;

		return result;
	}

	inline mat4 identity() {
		return mat4(1.f);
	}

	inline mat4 perspective(real left, real right, real bottom, real top, real near, real far) {

		mat4 result = mat4();

		result.data[0][0]  = 2 * near / (right - left);
		result.data[0][2]  = (right + left) / (right - left);
		result.data[1][1]  = 2 * near / (top - bottom);
		result.data[1][2]  = (top + bottom) / (top - bottom);
		result.data[2][2] = -(far + near) / (far - near);
		result.data[2][3] = -(2 * far * near) / (far - near);
		result.data[3][2] = -1;
		result.data[3][3] = 0;

		return result;
	}

	inline mat4 perspective(real fov, real aspect, real near, real far) {

		real height = near * tan(radians(fov / 2.f));
		real width = height * aspect;

		return perspective(-width, width, -height, height, near, far);
	}

	inline mat4 ortho(real left, real right, real bottom, real top, real near, real far) {

		mat4 result = mat4();

		result.data[0][0]  = 2 / (right - left);
		result.data[0][3]  = -(right + left) / (right - left);
		result.data[1][1]  = 2 / (top - bottom);
		result.data[1][3]  = -(top + bottom) / (top - bottom);
		result.data[2][2] = -2 / (far - near);
		result.data[2][3] = -(far + near) / (far - near);

		return result;
	}


}

#endif
