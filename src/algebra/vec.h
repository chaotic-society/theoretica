#ifndef UROBORO_VECTOR_H
#define UROBORO_VECTOR_H
#include <array>

#include "../real_analysis.h"
#include "../vec_buff.h"

namespace uroboro {

	template<unsigned int N, typename T = real>
	class vec {

		public:

		static constexpr unsigned int size = N;

		T data[N];


		vec() = default;

		// Copy constructor
		vec(const vec<N, T>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.data[i];
			}
		}

		// Copy from other
		vec<N, T>& operator=(const vec<N, T>& other) {
			for (int i = 0; i < N; ++i) {
				data[i] = other.data[i];
			}
			return *this;
		}

		// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<T> l) {

			if(l.size() != N)
				// throw ...
				return;

			std::copy(l.begin(), l.end(), &data[0]);
		}

		~vec() = default;

		// Operators

		// Vector sum (v + w = (v.x + w.x, ...))
		inline vec<N, T> operator+(const vec<N, T>& other) const {
			vec<N, T> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] + other.data[i];
			}

			return result;
		}

		inline vec<N, T> operator-(const vec<N, T>& other) const {
			vec<N, T> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] - other.data[i];
			}

			return result;
		}

		// Scalar multiplication (av = (v.x * a, ...))
		inline vec<N, T> operator*(T scalar) const {
			vec<N, T> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = scalar * data[i];
			}

			return result;
		}

		// Scalar division (v / a = (v.x / a, ...))
		inline vec<N, T> operator/(T scalar) const {
			vec<N, T> result;

			for (int i = 0; i < N; ++i) {
				result.data[i] = data[i] / scalar;
			}

			return result;
		}

		// Dot product between vectors (v * w = v.x * w.x + ...)
		inline T dot(const vec<N, T>& other) const {

			T result = 0;

			for (int i = 0; i < N; ++i) {
				result += data[i] * other.data[i];
			}

			return result;
		}

		// Dot product between vectors (v * w = v.x * w.x + ...)
		inline T operator*(const vec<N, T>& other) const {
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


		inline void operator+=(const vec<N, T>& other) {

			for (int i = 0; i < N; ++i)
				data[i] += other.data[i];
		}

		inline void operator-=(const vec<N, T>& other) {

			for (int i = 0; i < N; ++i)
				data[i] -= other.data[i];
		}

		inline void operator*=(T scalar) {

			for (int i = 0; i < N; ++i)
				data[i] *= scalar;
		}

		inline void operator/=(T scalar) {

			for (int i = 0; i < N; ++i)
				data[i] /= scalar;
		}


		// Magnitude of vector (sqrt(v * v))
		inline T magnitude() const {

			T m = 0;
			for (int i = 0; i < N; ++i)
				m += data[i] * data[i];

			return sqrt(m);
		}

		// Alias for magnitude()
		inline T lenght() {
			return magnitude();
		}

		// Square magnitude of vector (v * v)
		inline T square_magnitude() const {
			T m = 0;
			for (int i = 0; i < N; ++i) {
				m += data[i] * data[i];
			}
			return m;
		}

		// Alias for square_magnitude()
		inline real square_lenght() const {
			return square_magnitude();
		}

		// Access i-th component
		inline T& operator[](int i) {
			return data[i];
		}

		// Access i-th element
		inline T& at(int i) {
			return data[i];
		}

		// Getters and setters
		inline T get(int i) const {
			return data[i];
		}

		// Set the i-th element
		inline void set(int i, T x) {
			data[i] = x;
		}

		// Vector normalization (v / |v|)
		inline void normalize() {
			real m = magnitude();

			if(m == 0)
				// throw ...
				return;

			for (int i = 0; i < N; ++i) {
				data[i] /= m;
			}
		}

		// Return the normalized vector (v / |v|)
		inline vec<N, T> normalized() const {
			vec<N, T> result = vec<N, T>();
			real m = magnitude();

			if(m == 0)
				// throw ...
				return result;

			for (int i = 0; i < N; ++i)
				result[i] = data[i] / m;

			return result;
		}

		// Check whether all elements of both vectors are equal
		inline bool operator==(const vec<N, T>& other) const {
			for (int i = 0; i < N; ++i) {
				if(data[i] != other[i])
					return false;
			}

			return true;
		}

		// Convert a vec<N, T> to a vec_buff
		inline vec_buff to_vec_buff() {
			
			vec_buff res;
			for (int i = 0; i < N; ++i)
				res.push_back(static_cast<real>(data[i]));

			return res;
		}

	};

	// Common vector types
	using vec2 = vec<2, real>;
	using vec3 = vec<3, real>;
	using vec4 = vec<4, real>;


	// Compute the (Euclidian) distance between two values
	real distance(real v1, real v2) {
		return abs(v1 - v2);
	}


	// Compute the distance between two points
	template<unsigned int N, typename T>
	T distance(vec<N, T> v1, vec<N, T> v2) {
		return (v1 - v2).lenght();
	}


	// template<unsigned int N>
	// real lp_norm(vec<N, T> v1, vec<N, T> v2, unsigned int order) {

	// 	real sum = 0;

	// 	for (unsigned int i = 0; i < N; ++i) {
	// 		sum += uroboro::pow(v1[i] - v2[i], order);
	// 	}

	// 	return uroboro::powf(sum, 1.0 / (real) order);
	// }


	// Compute the dot product of two vectors
	template<unsigned int N, typename T>
	real dot(vec<N, T> v1, vec<N, T> v2) {
		return v1.dot(v2);
	}


	// Compute the cross product of two vectors
	vec3 cross(const vec3& v1, const vec3& v2) {
		vec3 res;

		res.data[0] = v1.data[1] * v2.data[2] - v1.data[2] * v2.data[1];
		res.data[1] = v1.data[2] * v2.data[0] - v1.data[0] * v2.data[2];
		res.data[2] = v1.data[0] * v2.data[1] - v1.data[1] * v2.data[0];

		return res;
	}

}


#endif
