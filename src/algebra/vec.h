
/// 
/// @file vec.h Vector class and operations
/// 

#ifndef THEORETICA_VECTOR_H
#define THEORETICA_VECTOR_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/error.h"
#include "../core/real_analysis.h"
#include "./algebra.h"
#include <vector>


namespace theoretica {

	
	/// @class vec_iterator
	/// A sequential iterator for traversing vector-like containers.
	///
	/// @tparam Vector The type of vector being iterated over.
	/// @tparam ReturnType The type returned when dereferencing the iterator, defaults to a reference to the vector's element type.
	template<typename Vector, typename ReturnType = vector_element_t<Vector>&>
	class vec_iterator {

		private:

			/// A reference to the vector being iterated over
			Vector& vector;

			/// The current index
			size_t i;
		
		public:

			using iterator_category = std::forward_iterator_tag;
			using value_type = vector_element_t<Vector>;
			using pointer = value_type*;
			using reference = value_type&;

			/// Construct the iterator from a pointer to the
			/// elements and a starting index.
			vec_iterator(Vector& vector, size_t index) : vector(vector), i(index) {}

			/// Dereference the iterator
			/// to get the current element.
			ReturnType operator*() {
				return vector[i];
			}

			/// Move to the next element in the vector.
			vec_iterator& operator++() {
				++i;
				return *this;
			}

			/// Move to the previous element in the vector.
			// vec_iterator& operator--() {
			// 	--i;
			// 	return *this;
			// }


			/// Get the current index
			size_t index() {
				return i;
			}


			/// Comparison operators.
			bool operator==(const vec_iterator& other) const {
				return i == other.i;
			}

			bool operator!=(const vec_iterator& other) const {
				return !(*this == other);
			}
	};


	/// 
	/// @class vec
	/// A statically allocated N-dimensional vector
	/// with elements of the given type.
	/// 
	template<typename Type = real, unsigned int N = 0>
	class vec {
		
		private:
			Type elements[N];

		public:

		// Vector size template argument
		static constexpr size_t size_argument = N;

		/// Construct a vector with all
		/// elements equal to zero.
		vec() {
			algebra::vec_zeroes(*this);
		}


		/// Construct a vector with all elements
		/// equal to the given value.
		vec(Type val) {

			for (unsigned int i = 0; i < N; ++i)
				elements[i] = val;
		}


		/// Construct a vector with all elements
		/// equal to the given value, checking that
		/// the given size matches that of the vector type.
		vec(unsigned int size, Type val) {

			if (N != size) {
				TH_MATH_ERROR("vec::vec(size, val)", size, MathError::InvalidArgument);
				algebra::vec_error(*this);
				return;
			}

			for (unsigned int i = 0; i < N; ++i)
				elements[i] = val;
		}


		/// Copy constructor
		template<unsigned int M>
		vec(const vec<Type, M>& other) {

			if (N == other.size())
				algebra::vec_copy(*this, other);
			else
				algebra::vec_error(*this);
		}


		/// Construct a vector from its elements, provided they are more than two
		/// (to avoid conflict with other constructors).
		template<typename... Args>
		vec(Type x1, Type x2, Args... args) {

			static_assert(
				2 + sizeof...(args) == N,
				"Number of arguments must match vector size"
			);

			*this = {x1, x2, args...};
		}


		/// Copy from other
		template<typename Vector>
		vec<Type, N>& operator=(const Vector& other) {

			if (N == other.size())
				return algebra::vec_copy(*this, other);
			else
				return algebra::vec_error(*this);
		}


		/// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<Type> l) {

			if(l.size() != N) {
				TH_MATH_ERROR("vec::vec(initializer_list<Type>)", l.size(), MathError::InvalidArgument);
				algebra::vec_error(*this);
				return;
			}

			std::copy(l.begin(), l.end(), &elements[0]);
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
				result.elements[i] = scalar * elements[i];

			return result;
		}


		/// Scalar division (v / a = (v.x / a, ...))
		inline vec<Type, N> operator/(Type scalar) const {
			vec<Type, N> result;

			for (unsigned int i = 0; i < N; ++i)
				result.elements[i] = elements[i] / scalar;

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
			return algebra::cross(*this, other);
		}


		/// Sum a vector the the vector itself
		template<typename Vector>
		inline vec<Type, N>& operator+=(const Vector& other) {

			for (unsigned int i = 0; i < N; ++i)
				elements[i] += other.elements[i];
		
			return *this;
		}


		/// Subtract a vector the the vector itself
		template<typename Vector>
		inline vec<Type, N>& operator-=(const Vector& other) {

			for (unsigned int i = 0; i < N; ++i)
				elements[i] -= other.elements[i];
		
			return *this;
		}


		/// Multiply the vector itself by a scalar
		inline vec<Type, N>& operator*=(Type scalar) {

			for (unsigned int i = 0; i < N; ++i)
				elements[i] *= scalar;
		
			return *this;
		}


		/// Divide the vector itself by a scalar
		inline vec<Type, N>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("vec::operator/=", scalar, MathError::DivByZero);
				algebra::vec_error(*this);
				return *this;
			}

			for (unsigned int i = 0; i < N; ++i)
				elements[i] /= scalar;
		
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


		/// Access i-th component by reference.
		inline Type& operator[](unsigned int i) {
			return elements[i];
		}


		/// Get the i-th component by value.
		inline const Type& operator[](unsigned int i) const {
			return elements[i];
		}


		/// Access i-th element by reference, with bound checking.
		///
		/// If the given index is out of range, an std::out_of_range
		/// exception is thrown.
		inline Type& at(unsigned int i) {

			if (i >= N) {
				throw std::out_of_range(
					"The element index in vec::at() is out of bounds"
				);
			}

			return elements[i];
		}


		/// Get the i-th element by value, with bound checking.
		///
		/// If the given index is out of range, an std::out_of_range
		/// exception is thrown.
		inline Type at(unsigned int i) const {

			if (i >= N) {
				throw std::out_of_range(
					"The element index in vec::at() is out of bounds"
				);
			}

			return elements[i];
		}


		/// Sequential iterator for statically allocated vectors.
		using iterator = vec_iterator<vec<Type, N>, Type&>;

		/// Const sequential iterator for statically allocated vectors.
		using const_iterator = vec_iterator<const vec<Type, N>, const Type&>;


		/// Get an iterator to the first element
		/// of the vector.
		inline auto begin() {
			return iterator(*this, 0);
		}


		/// Get an iterator to one plus the last element
		/// of the vector.
		inline auto end() {
			return iterator(*this, size());
		}


		/// Get a const iterator to the first element of the vector.
		inline auto begin() const {
			return const_iterator(*this, 0);
		}


		/// Get a const iterator to one plus the
		/// last element of the vector.
		inline auto end() const {
			return const_iterator(*this, size());
		}


		/// Get a raw pointer to the elements of the vector.
		inline Type* data() {
			return elements;
		}


		/// Get a raw pointer to the elements of the vector.
		inline const Type* data() const {
			return elements;
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
				if(elements[i] != other[i])
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
				TH_MATH_ERROR("vec::resize", N, MathError::InvalidArgument);
			}
		}


		/// Compatibility function to allow for allocation
		/// or resizing of dynamic vectors. Since statically
		/// allocated vectors cannot change size, this function
		/// only checks whether the target size is the same
		/// as the vector's.
		///
		/// The value argument is ignored since the size cannot be changed.
		inline void resize(size_t n, const Type& value) const {
			
			if(N != n) {
				TH_MATH_ERROR("vec::resize", N, MathError::InvalidArgument);
			}
		}


		/// Returns an N-dimensional euclidean base unit vector
		/// with the i-th element set to 1.
		inline static vec<Type, N> euclidean_base(
			unsigned int i, unsigned int n = N) {

			if(i >= n) {
				TH_MATH_ERROR("vec::euclidean_base", i, MathError::InvalidArgument);
				return vec<Type, N>(Type(nan()));
			}

			vec<Type, N> e_i = vec<Type, N>(n, Type(0.0));
			e_i.resize(n);
			e_i[i] = 1;

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
				res << elements[i];
				if(i != N - 1)
					res << separator;
			}

			if(parenthesis)
				res << ")";

			return res.str();
		}


		/// Convert the vector to string representation.
		inline operator std::string() {
			return to_string();
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
	/// A dynamically allocated vector
	/// with elements of the given type.
	///
	template<typename Type>
	class vec<Type, 0> {

		// Container type for storage (alias for std::vector)
		template<typename T>
		using Container = std::vector<T>;

	private:
			Container<Type> elements;

	public:

		// Vector size template argument
		static constexpr size_t size_argument = 0;

		/// Construct an empty vector.
		vec() = default;

		/// Construct a vector with the given size
		/// and all elements equal to zero.
		vec(unsigned int n) {
			resize(n);
			algebra::vec_zeroes(*this);
		}

		/// Construct a vector with the given size
		/// and all elements equal to the given value
		vec(unsigned int n, Type a) {
			elements = std::vector<Type>(n, a);
		}

		/// Copy constructor
		template<unsigned int M>
		vec(const vec<Type, M>& other) {
			algebra::vec_copy(*this, other);
		}

		/// Construct a vector from its elements, provided they are more than two
		/// (to avoid conflict with other constructors).
		template<typename... Args>
		vec(Type x1, Type x2, Args... args) {

			elements = {x1, x2, args...};
		}

		
		/// Initialize from a list, e.g. {1, 2, 3}
		vec(std::initializer_list<Type> l) : elements(l) {}


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
				result.elements[i] = scalar * elements[i];

			return result;
		}


		/// Scalar division (v / a = (v.x / a, ...))
		inline vec<Type> operator/(Type scalar) const {

			vec<Type> result;
			result.resize(size());

			for (unsigned int i = 0; i < size(); ++i)
				result.elements[i] = elements[i] / scalar;

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
			return algebra::cross(*this, other);
		}


		/// Sum a vector to the vector itself
		template<typename Vector>
		inline vec<Type>& operator+=(const Vector& other) {

			// If the vector is uninitialized,
			// initialize it to be a zero vector
			// with the target size
			if(size() == 0)
				resize(other.size());

			if(size() != other.size()) {
				TH_MATH_ERROR("vec::operator+=", size(), MathError::InvalidArgument);
				return (*this = vec<Type>(max(size(), 1), nan()));
			}

			for (unsigned int i = 0; i < size(); ++i)
				elements[i] += other.elements[i];
		
			return *this;
		}


		/// Subtract a vector from the vector itself
		template<typename Vector>
		inline vec<Type>& operator-=(const Vector& other) {

			if(size() != other.size()) {
				TH_MATH_ERROR("vec::operator-=", size(), MathError::InvalidArgument);
				return (*this = vec<Type>(max(size(), 1), nan()));
			}

			for (unsigned int i = 0; i < size(); ++i)
				elements[i] -= other.elements[i];
		
			return *this;
		}


		/// Multiply the vector itself by a scalar
		inline vec<Type>& operator*=(Type scalar) {

			for (unsigned int i = 0; i < size(); ++i)
				elements[i] *= scalar;
		
			return *this;
		}


		/// Divide the vector itself by a scalar
		inline vec<Type>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("vec::operator/=", scalar, MathError::DivByZero);
				*this = vec<Type>(max(size(), 1), nan());
				return *this;
			}

			for (unsigned int i = 0; i < size(); ++i)
				elements[i] /= scalar;
		
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


		/// Access i-th component by reference.
		inline Type& operator[](unsigned int i) {
			return elements[i];
		}


		/// Get the i-th component by value.
		inline const Type& operator[](unsigned int i) const {
			return elements[i];
		}


		/// Access i-th element by reference, with bound checking.
		///
		/// If the given index is out of range, an std::MathError::OutOfRange
		/// exception is thrown.
		inline Type& at(unsigned int i) {
			return elements.at(i);
		}


		/// Get the i-th element by value, with bound checking.
		///
		/// If the given index is out of range, an std::MathError::OutOfRange
		/// exception is thrown.
		inline Type at(unsigned int i) const {
			return elements.at(i);
		}


		using iterator = typename Container<Type>::iterator;


		/// Get an iterator to the first element
		/// of the vector.
		inline auto begin() {
			return elements.begin();
		}


		/// Get a const iterator to the first element
		/// of the vector.
		inline auto begin() const {
			return elements.cbegin();
		}


		/// Get an iterator to one plus the last element
		/// of the vector.
		inline auto end() {
			return elements.end();
		}


		/// Get a const iterator to one plus the last element
		/// of the vector.
		inline auto end() const {
			return elements.cend();
		}


		/// Get a raw pointer to the elements of the vector.
		inline Type* data() {
			return elements.data();
		}


		/// Get a raw pointer to the elements of the vector.
		inline const Type* data() const {
			return elements.data();
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
				if(elements[i] != other[i])
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
			return elements.size();
		}


		/// Change the size of the vector
		inline void resize(size_t n) {
			elements.resize(n);
		}


		/// Change the size of the vector, filling new elements with the given value
		inline void resize(size_t n, const Type& value) {
			elements.resize(n, value);
		}


		/// Add a value at the end of the vector
		/// (only for dynamically allocated vectors).
		inline void append(const Type& x) {
			elements.push_back(x);
		}


		/// Add a value at the end of the vector
		/// (only for dynamically allocated vectors).
		inline void append(Type&& x) {
			elements.push_back(x);
		}


		/// Returns an euclidean base unit vector
		/// with the i-th element set to 1 and size n.
		inline static vec<Type> euclidean_base(
			unsigned int i, unsigned int n) {

			if(i >= n) {
				TH_MATH_ERROR("vec::euclidean_base", i, MathError::InvalidArgument);
				return vec<Type>(n, nan());
			}

			vec<Type> e_i = vec<Type>(n, Type(0.0));
			e_i.resize(n);
			e_i[i] = 1;

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
				res << elements[i];
				if(i != size() - 1)
					res << separator;
			}
			
			if(parenthesis)
				res << ")";

			return res.str();
		}


		/// Convert the vector to string representation.
		inline operator std::string() {
			return to_string();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const vec<Type>& obj) {
			return out << obj.to_string();
		}

#endif

	};

}


#endif
