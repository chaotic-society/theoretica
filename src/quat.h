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

			quat() {
				a = 0;
				v = vec3();
			}

			quat(real a, const vec3& v) {
				this->a = a;
				this->v = v;
			}

			quat(const quat& other) {
				a = other.a;
				v = other.v;
			}

			inline quat& operator=(const quat& other) {
				a = other.a;
				v = other.v;
				return *this;
			}

			quat(real a, real b, real c, real d) {
				this-> a = a;
				v.data[0] = b;
				v.data[1] = c;
				v.data[2] = d;
			}

			~quat() {}

			// Operators

			// Quaternion sum


	};

}

#endif
