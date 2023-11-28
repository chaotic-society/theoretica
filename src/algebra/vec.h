
/// 
/// @file vec.h Vector class and operations
/// 

#ifndef THEORETICA_VECTOR_H
#define THEORETICA_VECTOR_H

#include "../core/error.h"
#include "../core/real_analysis.h"
#include "../core/vec_buff.h"
#include "./algebra.h"


namespace theoretica {

	/// 
	/// @class vec
	/// An N-dimensional vector of Type elements.
	/// 
	template<typename Type = real, unsigned int N = 0>
	class vec {
		public:

		Type data[N];


		vec() {
			algebra::vec_zeroes(*this);
		}


		/// Initialize all elements to the same value
		vec(Type a) {
			for (unsigned int i = 0; i < N; ++i)
				data[i] = a;
		}


		/// Copy constructor
		vec(const vec<Type, N>& other) {
			algebra::vec_copy(*this, other);
		}


		/// Copy from other
		template<typename Vector>
		vec<Type, N>& operator=(const Vector& other) {
			return algebra::vec_copy(*this, other);
		}


		/// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<Type> l) {

			if(l.size() != N) {
				TH_MATH_ERROR("vec::vec(initializer_list<Type>)", l.size(),
					INVALID_ARGUMENT);
				// Set all elements to NaN
				*this = vec<Type, N>(nan());
				return;
			}

			std::copy(l.begin(), l.end(), &data[0]);
		}

		~vec() = default;


		/// Identity
		inline vec<Type, N> operator+() const {
			return *this;
		}


		/// Vector sum (v + w = (v.x + w.x, ...))
		inline vec<Type, N> operator+(const vec<Type, N>& other) const {
			
			vec<Type, N> result;
			algebra::vec_sum(result, *this, other);
			return result;
		}


		/// Opposite vector
		inline vec<Type, N> operator-() const {
			return *this * (Type) -1;
		}


		/// Vector subtraction
		inline vec<Type, N> operator-(const vec<Type, N>& other) const {
			
			vec<Type, N> result;
			algebra::vec_diff(result, *this, other);
			return result;
		}


		/// Scalar multiplication (av = (v.x * a, ...))
		inline vec<Type, N> operator*(Type scalar) const {
			vec<Type, N> result;

			for (unsigned int i = 0; i < N; ++i)
				result.data[i] = scalar * data[i];

			return result;
		}


		/// Scalar division (v / a = (v.x / a, ...))
		inline vec<Type, N> operator/(Type scalar) const {
			vec<Type, N> result;

			for (unsigned int i = 0; i < N; ++i)
				result.data[i] = data[i] / scalar;

			return result;
		}


		/// Dot product between vectors (v * w = v.x * w.x + ...)
		template<typename Vector>
		inline Type dot(const Vector& other) const {
			return algebra::dot(*this, other);
		}


		/// Dot product between vectors (v * w = v.x * w.x + ...)
		template<typename Vector>
		inline Type operator*(const Vector& other) const {
			return dot(other);
		}


		/// Cross product between vectors
		inline vec<Type, N> cross(const vec<Type, N>& other) const {
			static_assert(N == 3, "The vector must be three dimensional");
			return algebra::cross(*this, other);
		}


		/// Cross product between vectors
		template<typename Vector>
		inline vec<Type, N> cross(const Vector& other) const {
		
			if(other.size() != 3) {
				TH_MATH_ERROR("vec::cross", other.size(), INVALID_ARGUMENT);
				return vec<Type, N>((Type) nan());
			}

			return algebra::cross(*this, other);
		}


		/// Sum a vector the the vector itself
		template<typename Vector>
		inline vec<Type, N>& operator+=(const Vector& other) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] += other.data[i];
		
			return *this;
		}


		/// Subtract a vector the the vector itself
		template<typename Vector>
		inline vec<Type, N>& operator-=(const Vector& other) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] -= other.data[i];
		
			return *this;
		}


		/// Multiply the vector itself by a scalar
		inline vec<Type, N>& operator*=(Type scalar) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] *= scalar;
		
			return *this;
		}


		/// Divide the vector itself by a scalar
		inline vec<Type, N>& operator/=(Type scalar) {

			for (unsigned int i = 0; i < N; ++i)
				data[i] /= scalar;
		
			return *this;
		}


		/// Compute the norm of the vector (sqrt(v * v))
		inline Type norm() const {
			return algebra::norm(*this);
		}


		/// Compute the square norm of the vector (v * v)
		inline Type sqr_norm() const {
			return algebra::sqr_norm(*this);
		}


		/// Access i-th component
		inline Type& operator[](unsigned int i) {
			return data[i];
		}


		/// Get the i-th component
		inline const Type& operator[](unsigned int i) const {
			return data[i];
		}


		/// Access i-th element
		inline Type& at(unsigned int i) {
			return data[i];
		}


		/// Getters and setters
		inline Type get(unsigned int i) const {
			return data[i];
		}


		/// Set the i-th element
		inline void set(unsigned int i, Type x) {
			data[i] = x;
		}


		/// Vector normalization (v / |v|)
		inline void normalize() {
			algebra::make_normalized(*this);
		}


		/// Return the normalized vector (v / |v|)
		inline vec<Type, N> normalized() const {
			return algebra::normalize(*this);
		}


		/// Check whether all elements of both vectors are equal
		template<typename Vector>
		inline bool operator==(const Vector& other) const {

			if(size() != other.size())
				return false;

			for (unsigned int i = 0; i < N; ++i)
				if(data[i] != other[i])
					return false;

			return true;
		}


		/// Check whether all elements of both vectors are unequal
		template<typename Vector>
		inline bool operator!=(const Vector& other) const {
			return !(*this == other);
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
		inline static vec<Type, N> euclidean_base(unsigned int i) {

			vec<Type, N> e_i = vec<Type, N>(0);
			e_i.at(i) = 1;

			return e_i;
		}


		/// Friend operator to enable equations of the form
		/// (Type) * (vec)
		inline friend vec<Type, N> operator*(Type a, const vec<Type, N>& v) {
			return v * a;
		}


#ifndef THEORETICA_NO_PRINT

		/// Convert the vector to string representation
		inline std::string to_string(
			const std::string& separator = ", ",
			bool parenthesis = true) const {

			std::stringstream res;

			if(parenthesis)
				res << "(";
			
			for (unsigned int i = 0; i < N; ++i) {
				res << data[i];
				if(i != N - 1)
					res << separator;
			}

			if(parenthesis)
				res << ")";

			return res.str();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(
			std::ostream& out, const vec<Type, N>& obj) {
			return out << obj.to_string();
		}

#endif

	};


	/// 
	/// @class vec
	/// A vector of Type elements and variable size.
	/// 
	template<typename Type>
	class vec<Type, 0> {
		public:

		std::vector<Type> data;

		vec() {}


		/// Initialize to the given size
		vec(unsigned int n) {
			resize(n);
			algebra::vec_zeroes(*this);
		}


		/// Initialize all elements to the same value
		vec(Type a, unsigned int n) {
			data = std::vector<Type>(a, n);
		}


		/// Copy constructor
		template<unsigned int N>
		vec(const vec<Type, N>& other) {
			algebra::vec_copy(*this, other);
		}


		/// Copy from other
		template<typename Vector>
		vec<Type>& operator=(const Vector& other) {
			return algebra::vec_copy(*this, other);
		}


		/// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<Type> l) {

			if(l.size() != size()) {
				TH_MATH_ERROR("vec::vec(initializer_list<Type>)", l.size(),
					INVALID_ARGUMENT);
				// Set all elements to NaN
				*this = vec<Type>(nan(), size());
				return;
			}

			std::copy(l.begin(), l.end(), &data[0]);
		}

		~vec() = default;


		/// Identity
		inline vec<Type> operator+() const {
			return *this;
		}


		/// Vector sum (v + w = (v.x + w.x, ...))
		template<typename Vector>
		inline vec<Type> operator+(const Vector& other) const {
			
			vec<Type> result;
			result.resize(size());
			algebra::vec_sum(result, *this, other);
			return result;
		}


		/// Opposite vector
		inline vec<Type> operator-() const {
			return *this * (Type) -1;
		}


		/// Vector subtraction
		template<typename Vector>
		inline vec<Type> operator-(const Vector& other) const {
			
			vec<Type> result;
			result.resize(size());
			algebra::vec_diff(result, *this, other);
			return result;
		}


		/// Scalar multiplication (av = (v.x * a, ...))
		inline vec<Type> operator*(Type scalar) const {

			vec<Type> result;
			result.resize(size());

			for (unsigned int i = 0; i < size(); ++i)
				result.data[i] = scalar * data[i];

			return result;
		}


		/// Scalar division (v / a = (v.x / a, ...))
		inline vec<Type> operator/(Type scalar) const {

			vec<Type> result;
			result.resize(size());

			for (unsigned int i = 0; i < size(); ++i)
				result.data[i] = data[i] / scalar;

			return result;
		}


		/// Dot product between vectors (v * w = v.x * w.x + ...)
		template<typename Vector>
		inline Type dot(const Vector& other) const {
			return algebra::dot(*this, other);
		}


		/// Dot product between vectors (v * w = v.x * w.x + ...)
		template<typename Vector>
		inline Type operator*(const Vector& other) const {
			return dot(other);
		}


		/// Cross product between vectors
		template<typename Vector>
		inline vec<Type> cross(const Vector& other) const {
		
			if(other.size() != 3) {
				TH_MATH_ERROR("vec::cross", other.size(), INVALID_ARGUMENT);
				return vec<Type>((Type) nan(), 3);
			}

			return algebra::cross(*this, other);
		}


		/// Sum a vector to the vector itself
		template<typename Vector>
		inline vec<Type>& operator+=(const Vector& other) {

			if(size() != other.size()) {
				TH_MATH_ERROR("vec::operator+=", size(), INVALID_ARGUMENT);
				return vec<Type>(nan(), size());
			}

			for (unsigned int i = 0; i < size(); ++i)
				data[i] += other.data[i];
		
			return *this;
		}


		/// Subtract a vector the the vector itself
		template<typename Vector>
		inline vec<Type>& operator-=(const Vector& other) {

			if(size() != other.size()) {
				TH_MATH_ERROR("vec::operator-=", size(), INVALID_ARGUMENT);
				return vec<Type>(nan(), size());
			}

			for (unsigned int i = 0; i < size(); ++i)
				data[i] -= other.data[i];
		
			return *this;
		}


		/// Multiply the vector itself by a scalar
		inline vec<Type>& operator*=(Type scalar) {

			for (unsigned int i = 0; i < size(); ++i)
				data[i] *= scalar;
		
			return *this;
		}


		/// Divide the vector itself by a scalar
		inline vec<Type>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("vec::operator/=", scalar, DIV_BY_ZERO);
				return vec<Type>(nan(), size());
			}

			for (unsigned int i = 0; i < size(); ++i)
				data[i] /= scalar;
		
			return *this;
		}


		/// Compute the norm of the vector (sqrt(v * v))
		inline Type norm() const {
			return algebra::norm(*this);
		}


		/// Compute the square norm of the vector (v * v)
		inline Type sqr_norm() const {
			return algebra::sqr_norm(*this);
		}


		/// Access i-th component
		inline Type& operator[](unsigned int i) {
			return data[i];
		}


		/// Get the i-th component
		inline const Type& operator[](unsigned int i) const {
			return data[i];
		}


		/// Access i-th element
		inline Type& at(unsigned int i) {
			return data[i];
		}


		/// Getters and setters
		inline Type get(unsigned int i) const {
			return data[i];
		}


		/// Set the i-th element
		inline void set(unsigned int i, Type x) {
			data[i] = x;
		}


		/// Vector normalization (v / |v|)
		inline void normalize() {
			algebra::make_normalized(*this);
		}


		/// Return the normalized vector (v / |v|)
		inline vec<Type> normalized() const {
			return algebra::normalize(*this);
		}


		/// Check whether all elements of both vectors are equal
		template<typename Vector>
		inline bool operator==(const Vector& other) const {

			if(size() != other.size())
				return false;

			for (unsigned int i = 0; i < size(); ++i)
				if(data[i] != other[i])
					return false;

			return true;
		}


		/// Check whether all elements of both vectors are unequal
		template<typename Vector>
		inline bool operator!=(const Vector& other) const {
			return !(*this == other);
		}


		/// Returns the size of the vector
		inline TH_CONSTEXPR unsigned int size() const {
			return data.size();
		}


		/// Compatibility function to allow for allocation
		/// or resizing of dynamic vectors. Since statically
		/// allocated vectors cannot change size, this function
		/// only checks whether the target size is the same
		/// as the vector's.
		inline void resize(size_t n) {
			data.resize(n);
		}


		/// Returns an euclidean base unit vector
		/// with the i-th element set to 1 and size n.
		inline static vec<Type> euclidean_base(
			unsigned int i, unsigned int n) {

			vec<Type> e_i = vec<Type>(0, n);
			e_i.at(i) = 1;

			return e_i;
		}


		/// Friend operator to enable equations of the form
		/// (Type) * (vec)
		inline friend vec<Type> operator*(Type a, const vec<Type>& v) {
			return v * a;
		}


#ifndef THEORETICA_NO_PRINT

		/// Convert the vector to string representation
		inline std::string to_string(
			const std::string& separator = ", ",
			bool parenthesis = true) const {

			std::stringstream res;

			if(parenthesis)
				res << "(";
			
			for (unsigned int i = 0; i < size(); ++i) {
				res << data[i];
				if(i != size() - 1)
					res << separator;
			}
			
			if(parenthesis)
				res << ")";

			return res.str();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const vec<Type>& obj) {
			return out << obj.to_string();
		}

#endif

	};

}


#endif
