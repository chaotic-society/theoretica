#ifndef UROBORO_VECTOR_H
#define UROBORO_VECTOR_H

#include <array>
#include "./common.h"

namespace uroboro {

	template<unsigned int N>
	class vec {
		public:

		const unsigned int size = N;

		real data[N];

		inline vec() {
			for (unsigned int i = 0; i < N; ++i) {
				data[i] = 0;
			}
		}

		inline vec(const vec<N>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.data[i];
			}
		}

		inline vec(const std::array<real, N>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.at(i);
			}
		}

		inline ~vec() {}

		inline vec<N> operator+(const vec<N>& other) {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] + other.data[i];
			}

			return result;
		}

		inline vec<N> operator-(const vec<N>& other) {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] - other.data[i];
			}

			return result;
		}

		inline vec<N> operator*(real scalar) {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = scalar * data[i];
			}

			return result;
		}

		inline real dot(const vec<N>& other) {

			real result = 0;

			for (int i = 0; i < N; ++i) {
				result += data[i] * other.data[i];
			}

			return result;
		}

		inline real magnitude() {
			real m = 0;
			for (int i = 0; i < N; ++i) {
				m += data[i] * data[i];
			}
			return sqrt(m);
		}

		inline real lenght() {
			return magnitude();
		}

		inline real square_magnitude() {
			real m = 0;
			for (int i = 0; i < N; ++i) {
				m += data[i] * data[i];
			}
			return m;
		}

		inline real square_lenght() {
			return square_magnitude();
		}

		inline real& operator[](int i) {
			return data[i];
		}

		inline void normalize() {
			real m = magnitude();

			if(m == 0)
				return;

			for (int i = 0; i < N; ++i) {
				data[i] /= m;
			}
		}

		inline vec<N> normalized() {
			vec<N> result = vec<N>();
			real m = magnitude();

			if(m == 0)
				return result;

			for (int i = 0; i < N; ++i) {
				result[i] = data[i] / m;
			}

			return result;
		}

		inline bool operator==(const vec<N>& other) {
			for (int i = 0; i < N; ++i) {
				if(data[i] != other[i])
					return false;
			}

			return true;
		}

	};

}


#endif
