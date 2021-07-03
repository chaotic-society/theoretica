#ifndef UROBORO_VECTOR_H
#define UROBORO_VECTOR_H
#include <array>

#include "./common.h"
#include "./vec_buff.h"

namespace uroboro {

	template<unsigned int N>
	class vec {

		public:

		static constexpr unsigned int size = N;

		real data[N];


		vec() = default;

		// Copy constructor
		vec(const vec<N>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.data[i];
			}
		}

		// Copy from other
		vec<N>& operator=(const vec<N>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.data[i];
			}
			return *this;
		}

		// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<real> l) {

			if(l.size() != N)
				// throw ...
				return;

			std::copy(l.begin(), l.end(), &data[0]);
		}

		~vec() = default;

		// Operators

		// Vector sum (v + w = (v.x + w.x, ...))
		inline vec<N> operator+(const vec<N>& other) const {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] + other.data[i];
			}

			return result;
		}

		inline vec<N> operator-(const vec<N>& other) const {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] - other.data[i];
			}

			return result;
		}

		// Scalar multiplication (av = (v.x * a, ...))
		inline vec<N> operator*(real scalar) const {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = scalar * data[i];
			}

			return result;
		}

		inline vec<N> operator/(real scalar) const {
			vec<N> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] / scalar;
			}

			return result;
		}

		// Dot product between vectors (v * w = v.x * w.x + ...)
		inline real dot(const vec<N>& other) const {

			real result = 0;

			for (int i = 0; i < N; ++i) {
				result += data[i] * other.data[i];
			}

			return result;
		}

		inline real operator*(const vec<N>& other) const {
			return dot(other);
		}

		// Cross product between vectors
		inline vec<3> cross(const vec<3>& other) const {

			if(N != 3)
				//throw ...
				return vec<3>();

			vec<3> res;

			res.data[0] = data[1] * other.data[2] - data[2] * other.data[1];
			res.data[1] = data[2] * other.data[0] - data[0] * other.data[2];
			res.data[2] = data[0] * other.data[1] - data[1] * other.data[0];

			return res;
		}

		// wedge product (generalized cross product) ...


		inline void operator+=(const vec<N>& other) {

			for (int i = 0; i < N; ++i)
				data[i] += other.data[i];
		}

		inline void operator-=(const vec<N>& other) {

			for (int i = 0; i < N; ++i)
				data[i] -= other.data[i];
		}

		inline void operator*=(real scalar) {

			for (int i = 0; i < N; ++i)
				data[i] *= scalar;
		}

		inline void operator/=(real scalar) {

			for (int i = 0; i < N; ++i)
				data[i] /= scalar;
		}


		// Magnitude of vector (sqrt(v * v))
		inline real magnitude() const {

			real m = 0;
			for (int i = 0; i < N; ++i)
				m += data[i] * data[i];

			return sqrt(m);
		}

		inline real lenght() {
			return magnitude();
		}

		// Square magnitude of vector (v * v)
		inline real square_magnitude() const {
			real m = 0;
			for (int i = 0; i < N; ++i) {
				m += data[i] * data[i];
			}
			return m;
		}

		inline real square_lenght() const {
			return square_magnitude();
		}

		// Access i-th component
		inline real& operator[](int i) {
			return data[i];
		}

		inline real& at(int i) {
			return data[i];
		}

		// Getters and setters
		inline real get(int i) const {
			return data[i];
		}

		inline void set(int i, real x) {
			data[i] = x;
		}

		// Normalization (v / |v|)
		inline void normalize() {
			real m = magnitude();

			if(m == 0)
				// throw ...
				return;

			for (int i = 0; i < N; ++i) {
				data[i] /= m;
			}
		}

		inline vec<N> normalized() const {
			vec<N> result = vec<N>();
			real m = magnitude();

			if(m == 0)
				// throw ...
				return result;

			for (int i = 0; i < N; ++i)
				result[i] = data[i] / m;

			return result;
		}

		inline bool operator==(const vec<N>& other) const {
			for (int i = 0; i < N; ++i) {
				if(data[i] != other[i])
					return false;
			}

			return true;
		}

		// Convert a vec<N> to a vec_buff
		inline vec_buff to_vec_buff() {
			vec_buff res;
			for (int i = 0; i < N; ++i)
				res.push_back(data[i]);
			return res;
		}

	};

	// Common vector types
	using vec2 = vec<2>;
	using vec3 = vec<3>;
	using vec4 = vec<4>;

}


#endif
