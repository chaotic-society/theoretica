#ifndef _MAT_H
#define _MAT_H
#include "vec.h"

namespace uroboro {

	class mat4 {

		public:
			real data[4][4] = { {0.f}, {0.f}, {0.f}, {0.f} }; //4x4

			inline mat4() {}

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
				data[1][0] = 0;
				data[2][0] = 0;
				data[3][0] = 0;
				data[0][1] = 0;
				data[1][1] = 1;
				data[2][1] = 0;
				data[3][1] = 0;
				data[0][2] = 0;
				data[1][2] = 0;
				data[2][2] = 1;
				data[3][2] = 0;
				data[0][3] = 0;
				data[1][3] = 0;
				data[2][3] = 0;
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
				// Replace with fast memcpy
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
				return vec4((data[0][0] * vector.x) + (data[1][0] * vector.y) + (data[2][0] * vector.z) + data[3][0],
							(data[0][1] * vector.x) + (data[1][1] * vector.y) + (data[2][1] * vector.z) + data[3][1],
							(data[0][2] * vector.x) + (data[1][2] * vector.y) + (data[2][2] * vector.z) + data[3][2]);
			}

			inline void translate(vec4 translation) {
				data[3][0] += translation.x;
				data[3][1] += translation.y;
				data[3][2] += translation.z;
			}

			inline void rotate(real radians, vec4 rotation) {
			}

			//inline void rotate(quat rotation) {}

			inline void scale(vec4 scale) {
				data[0][0] *= scale.x;
				data[1][1] *= scale.y;
				data[2][2] *= scale.z;
			}

			inline void operator=(const mat4 &other) {
				// Replace with fast memcpy
				data[0][0] = other.data[0][0];
				data[1][0] = other.data[1][0];
				data[2][0] = other.data[2][0];
				data[3][0] = other.data[3][0];
				data[0][1] = other.data[0][1];
				data[1][1] = other.data[1][1];
				data[2][1] = other.data[2][1];
				data[3][1] = other.data[3][1];
				data[0][2] = other.data[0][2];
				data[1][2] = other.data[1][2];
				data[2][2] = other.data[2][2];
				data[3][2] = other.data[3][2];
				data[0][3] = other.data[0][3];
				data[1][3] = other.data[1][3];
				data[2][3] = other.data[2][3];
				data[3][3] = other.data[3][3];
			}

			//inline void operator=(const quat &quaternion) {}

			inline vec4 operator*(vec4 vector) {
				return vec4((data[0][0] * vector.x) + (data[1][0] * vector.y) + (data[2][0] * vector.z) + data[3][0],
							(data[0][1] * vector.x) + (data[1][1] * vector.y) + (data[2][1] * vector.z) + data[3][1],
							(data[0][2] * vector.x) + (data[1][2] * vector.y) + (data[2][2] * vector.z) + data[3][2]);
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

			inline real& operator[](unsigned int i) {
				return data[i][i];
			}

	};

}

#endif
