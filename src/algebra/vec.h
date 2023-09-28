
/// 
/// @file vec.h Vector class and operations
/// 

#ifndef THEORETICA_VECTOR_H
#define THEORETICA_VECTOR_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include <array>

#include "../core/error.h"
#include "../core/real_analysis.h"
#include "../core/vec_buff.h"
#include "./algebra.h"


namespace theoretica {

	/// 
	/// @class vec
	/// An N-dimensional vector of T type elements.
	/// 
	template<unsigned int N, typename T = real>
	class vec {
		public:

		static_assert(N > 0, "N cannot be zero.");
		T data[N];

		vec() = default;

		/// Initialize all elements to the same value
		vec(T a) {
			for (unsigned int i = 0; i < N; ++i) {
				data[i] = a;
			}
		}

		/// Copy constructor
		vec(const vec<N, T>& other) {
			algebra::vec_copy(*this, other);
		}

		/// Copy from other
		vec<N, T>& operator=(const vec<N, T>& other) {
			return algebra::vec_copy(*this, other);
		}

		/// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<T> l) {

			if(l.size() != N) {
				TH_MATH_ERROR("vec::vec(initializer_list<T>)", l.size(),
					INVALID_ARGUMENT);
				// Set all elements to NaN
				*this = vec<N, T>(nan());
				return;
			}

			std::copy(l.begin(), l.end(), &data[0]);
		}

		~vec() = default;


		// Operators

		/// Identity
		inline vec<N, T> operator+() const {
			return *this;
		}

		/// Vector sum (v + w = (v.x + w.x, ...))
		inline vec<N, T> operator+(const vec<N, T>& other) const {
			
			vec<N, T> result;
			algebra::vec_sum(result, *this, other);
			return result;
		}

		/// Opposite vector
		inline vec<N, T> operator-() const {
			return *this * -1;
		}

		inline vec<N, T> operator-(const vec<N, T>& other) const {
			
			vec<N, T> result;
			algebra::vec_diff(result, *this, other);
			return result;
		}

		/// Scalar multiplication (av = (v.x * a, ...))
		inline vec<N, T> operator*(T scalar) const {
			vec<N, T> result;

			for (unsigned int i = 0; i < N; ++i) {
				result.data[i] = scalar * data[i];
			}

			return result;
		}

		/// Scalar division (v / a = (v.x / a, ...))
		inline vec<N, T> operator/(T scalar) const {
			vec<N, T> result;

			for (unsigned int i = 0; i < N; ++i) {
				result.data[i] = data[i] / scalar;
			}

			return result;
		}


		/// Dot product between vectors (v * w = v.x * w.x + ...)
		inline T dot(const vec<N, T>& other) const {
			return algebra::dot(*this, other);
		}

		/// Dot product between vectors (v * w = v.x * w.x + ...)
		inline T operator*(const vec<N, T>& other) const {
			return dot(other);
		}


		/// Cross product between vectors
		inline vec<N> cross(const vec<N>& other) const {
			static_assert(N == 3, "The vector must be three dimensional");
			return algebra::cross(*this, other);
		}


		inline void operator+=(const vec<N, T>& other) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] += other.data[i];
		}

		inline void operator-=(const vec<N, T>& other) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] -= other.data[i];
		}

		inline void operator*=(T scalar) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] *= scalar;
		}

		inline void operator/=(T scalar) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] /= scalar;
		}


		/// Magnitude of vector (sqrt(v * v))
		inline T magnitude() const {
			return algebra::norm(*this);
		}

		/// Alias for magnitude()
		inline T length() const {
			return magnitude();
		}

		/// Square magnitude of vector (v * v)
		inline T square_magnitude() const {

			T m = 0;
			for (unsigned int i = 0; i < N; ++i)
				m += data[i] * conjugate(data[i]);

			return m;
		}


		/// Alias for square_magnitude()
		inline real square_length() const {
			return square_magnitude();
		}

		/// Access i-th component
		inline T& operator[](unsigned int i) {
			return data[i];
		}

		/// Access i-th element
		inline T& iat(unsigned int i) {
			return data[i];
		}

		/// Access i-th element
		inline T& at(unsigned int i) {
			return data[i];
		}

		/// Getters and setters
		inline T iget(unsigned int i) const {
			return data[i];
		}

		/// Getters and setters
		inline T get(unsigned int i) const {
			return data[i];
		}

		/// Set the i-th element
		inline void set(unsigned int i, T x) {
			data[i] = x;
		}


		/// Vector normalization (v / |v|)
		inline void normalize() {
			algebra::make_normalized(*this);
		}

		/// Return the normalized vector (v / |v|)
		inline vec<N, T> normalized() const {
			return algebra::normalize(*this);
		}


		/// Check whether all elements of both vectors are equal
		inline bool operator==(const vec<N, T>& other) const {
			for (unsigned int i = 0; i < N; ++i) {
				if(data[i] != other[i])
					return false;
			}

			return true;
		}


		/// Convert a vec<N, T> to a vec_buff
		inline vec_buff to_vec_buff() {
			
			vec_buff res;
			res.reserve(N);
			for (unsigned int i = 0; i < N; ++i)
				res.push_back(static_cast<real>(data[i]));

			return res;
		}


		/// Returns the size of the vector (N)
		inline TH_CONSTEXPR unsigned int size() const {
			return N;
		}


		/// Compatibility function to allow for allocation
		/// or resizing of dynamic vectors. Since statically
		/// allocated vectors cannot change size, this function
		/// only checks whether the target size is the same
		/// as the vector's.
		inline void resize(size_t n) const {
			
			if(N != n) {
				TH_MATH_ERROR("vec::resize", N, INVALID_ARGUMENT);
			}
		}


		/// Returns an N-dimensional euclidean base unit vector
		/// with the i-th element set to 1.
		inline static vec<N> euclidean_base(unsigned int i) {

			vec<N> e_i = vec<N>(0);
			e_i.at(i) = 1;

			return e_i;
		}

		/// Friend operator to enable equations of the form
		/// (T) * (vec)
		inline friend vec<N> operator*(T a, const vec<N, T>& v) {
			return v * a;
		}


#ifndef THEORETICA_NO_PRINT

		/// Convert the vector to string representation
		inline std::string to_string(const std::string& separator = ", ") const {

			std::stringstream res;

			res << "(";
			for (unsigned int i = 0; i < N; ++i) {
				res << data[i];
				if(i != N - 1)
					res << separator;
			}
			res << ")";

			return res.str();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const vec<N, T>& obj) {
			return out << obj.to_string();
		}

#endif

	};

}


#endif
