
///
/// @file mat.h Matrix class and operations
///

#ifndef THEORETICA_MATRIX_H
#define THEORETICA_MATRIX_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include <array>
#include <vector>

#include "../core/error.h"
#include "../core/constants.h"
#include "../core/real_analysis.h"
#include "./algebra.h"
#include "./transform.h"
#include "./vec.h"


namespace theoretica {


	/// @class mat_iterator
	/// A sequential iterator for matrices.
	///
	/// A const iterator may be constructed by specifying
	/// both typenames Matrix and ReturnType as const.
	template<typename Matrix, typename ReturnType = matrix_element_t<Matrix>&>
	class mat_iterator {

		private:

			/// Reference to the matrix being iterated over
			Matrix& matrix;

			/// Index of current row
			size_t row;

			/// Index of current column
			size_t col;

		public:

			using iterator_category = std::forward_iterator_tag;
			using value_type = matrix_element_t<Matrix>;
			using pointer = value_type*;
			using reference = value_type&;

        
			/// Constructs an iterator for a matrix, optionally starting at a specified row and column.
			/// @param matrix Reference to the matrix to iterate over.
			/// @param row Initial row index for the iterator (default is 0).
			/// @param col Initial column index for the iterator (default is 0).
			///
			/// Constructs an iterator that points to the specified initial position
			/// within the matrix. If no row and column are specified, the iterator
			/// starts at the matrix's top-left corner (0, 0).
			mat_iterator(
				Matrix& matrix,
				size_t row = 0,
				size_t col = 0)
			: matrix(matrix), row(row), col(col) {}
        	
			
			/// Dereferences the iterator to access the current matrix element.
			/// @return A reference (or const reference) to the matrix element at the current iterator position.
			///
			/// Provides direct access to the element at the iterator's current position.
			ReturnType operator*() {
				return matrix(row, col);
			}

        
			/// Advances the iterator to the next element in row-major order.
			/// @return A reference to the updated iterator.
			///
			/// Moves the iterator one element forward within the matrix. When the end
			/// of a row is reached, it continues to the first element of the next row.
			mat_iterator& operator++() {

				col++;

				if(col == matrix.cols()) {
					col = 0;
					row++;
				}

				return *this;
			}


			/// Retrieves the current row index of the iterator.
			/// @return The row index as a `size_t`.
			size_t row_index() {
				return row;
			}


			/// Retrieves the current column index of the iterator.
			/// @return The column index as a `size_t`.
			size_t col_index() {
				return col;
			}


			/// Equality operator to compare two iterators.
			/// @param other Another iterator to compare with.
			/// @return `true` if the iterators point to the same matrix element; otherwise, `false`.
			bool operator==(const mat_iterator& other) const {
				return (row == other.row) &&
					(col == other.col);
			}


			/// Inequality operator to compare two iterators.
			/// @param other Another iterator to compare with.
			/// @return `true` if the iterators do not point to the same matrix element; otherwise, `false`.
			bool operator!=(const mat_iterator& other) const {
				return !(*this == other);
			}
	};


	/// @class mat
	/// A generic matrix with a fixed number of rows and columns.
	///
	/// @param Type The type of the elements
	/// @param N The number of rows
	/// @param K The number of columns
	template<typename Type = real, unsigned int N = 0, unsigned int K = 0>
	class mat {
		public:

#ifdef THEORETICA_ROW_FIRST
		Type data[N][K];
#else
		Type data[K][N];
#endif


		/// Default constructor.
		///
		/// Initializes the matrix with all elements set to zero.
		mat() {
			algebra::mat_zeroes(*this);
		}


		/// Copy constructor.
		/// @tparam Matrix The type of the matrix to copy from.
		/// @param m The matrix to copy.
		///
		/// Copies all elements from the given matrix `m` into this matrix.
		template<typename Matrix>
		mat(const Matrix& m) {
			algebra::mat_copy(*this, m);
		}


		/// Constructs a matrix from an initializer list.
		/// @tparam T The type of elements in the initializer list (default is `Type`).
		/// @param rows An initializer list representing the rows of the matrix.
		///
		/// Initializes the matrix with values from the initializer list. Ensures that
		/// the initializer list matches the size of the matrix (i.e., `N` rows and `K` columns).
		template<typename T = Type>
		inline mat(const std::initializer_list<std::initializer_list<T>>& rows) {

			if(rows.size() != N) {
				TH_MATH_ERROR("mat::mat", rows.size(), INVALID_ARGUMENT);
				algebra::mat_error(*this);
				return;
			}

			unsigned int i = 0;
			unsigned int j = 0;

			for (const auto& row : rows) {

				if (row.size() != K) {
					TH_MATH_ERROR("mat::mat", rows.size(), INVALID_ARGUMENT);
					algebra::mat_error(*this);
					return;
				}

				for (const auto& x : row) {
					
					get(i, j) = x;
					j = (j + 1) % K;
				}

				i++;
			}
		}


		/// Constructor that initializes a diagonal matrix with equal entries on the diagonal.
		/// The size parameters are required for compatibility with mat<T, 0>
		/// and may be used for additional error prevention.
		///
		/// @param diagonal The value for the diagonal entries.
		/// @param n Number of rows.
		/// @param k Number of columns.
		mat(Type diagonal, unsigned int n = 0, unsigned int k = 0) {
			
			if (n && k)
				resize(n, k);

			algebra::mat_zeroes(*this);
			const unsigned int m = min(n, k);

			for (unsigned int i = 0; i < m; ++i)
				get(i, i) = diagonal;
		}


		/// Assignment operator to copy from another matrix.
		/// @tparam Matrix The type of the matrix to copy from.
		/// @param other The matrix to copy.
		/// @return Reference to this matrix after assignment.
		template<typename Matrix>
		inline mat<Type, N, K>& operator=(const Matrix& other) {
			return algebra::mat_copy(*this, other);
		}


		/// Matrix addition.
		/// @tparam Matrix The type of the matrix to add.
		/// @param other The matrix to add.
		/// @return A new matrix that is the result of adding `other` to this matrix.
		template<typename Matrix>
		inline mat<Type, N, K> operator+(const Matrix& other) const {
			mat<Type, N, K> res;
			return algebra::mat_sum(res, *this, other);
		}


		/// Matrix subtraction.
		/// @tparam Matrix The type of the matrix to subtract.
		/// @param other The matrix to subtract.
		/// @return A new matrix that is the result of subtracting `other` from this matrix.
		template<typename Matrix>
		inline mat<Type, N, K> operator-(const Matrix& other) const {
			mat<Type, N, K> res;
			return algebra::mat_diff(res, *this, other);
		}


		/// Scalar multiplication.
		/// @param scalar The scalar value to multiply.
		/// @return A new matrix that is the result of multiplying this matrix by `scalar`.
		inline mat<Type, N, K> operator*(Type scalar) const {
			mat<Type, N, K> res;
			return algebra::mat_scalmul(res, scalar, *this);
		}


		/// Friend operator for scalar multiplication (T * mat).
		/// @param a The scalar multiplier.
		/// @param B The matrix to be multiplied.
		/// @return A new matrix that is the result of multiplying `B` by `a`.
		inline friend mat<Type, N, K> operator*(Type a, const mat<Type, N, K>& B) {
			return B * a;
		}


		/// Friend operator for vector-matrix multiplication.
		/// @tparam VecType The type of vector elements.
		/// @tparam M The number of elements in the vector.
		/// @param a The vector to multiply.
		/// @param B The matrix to be multiplied by.
		/// @return The resulting vector after multiplication.
		template<typename VecType, unsigned int M>
		inline friend vec<VecType, K> operator*(
			const vec<VecType, M>& a, const mat<Type, N, K>& B) {
			return B * a;
		}


		/// Scalar division.
		/// @param scalar The scalar divisor.
		/// @return A new matrix that is the result of dividing this matrix by `scalar`.
		///
		/// If `scalar` is close to zero, an error is raised.
		inline mat<Type, N, K> operator/(Type scalar) const {

			mat<Type, N, K> res;

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transforms a vector by multiplying it with the matrix.
		/// @tparam Vector The type of the vector to transform.
		/// @param v The vector to transform.
		/// @return The transformed vector.
		///
		/// This function multiplies the given vector `v` by the matrix. It checks
		/// if the size of `v` matches the number of columns in the matrix.
		template<typename Vector>
		inline Vector transform(const Vector& v) const {

			if(v.size() != cols()) {
				TH_MATH_ERROR("mat::transform", v.size(), INVALID_ARGUMENT);
				Vector res;
				res.resize(cols());
				algebra::vec_error(res);
				return res;
			}

			return algebra::transform(*this, v);
		}


		/// Transforms a fixed-size vector by multiplying it with the matrix.
		/// @param v The vector to transform.
		/// @return The transformed vector as a new `vec<Type, N>`.
		inline vec<Type, N> transform(const vec<Type, K>& v) const {
			return algebra::transform(*this, v);
		}


   		/// Overloads the `*` operator to transform a fixed-size vector by the matrix.
   		/// @param v The vector to transform.
   		/// @return The transformed vector.
		inline vec<Type, N> operator*(const vec<Type, K>& v) const {
			return transform(v);
		}


		/// Matrix multiplication for matrices with different column counts.
		/// @tparam M The number of columns in matrix `B`.
		/// @param B The matrix to multiply with.
		/// @return A new matrix resulting from the multiplication.
		template<unsigned int M>
		inline mat<Type, N, M> mul(const mat<Type, K, M>& B) const {
			mat<Type, N, M> res;
			return algebra::mat_mul(res, *this, B);
		}


		/// Matrix multiplication for matrices with any type.
		/// @tparam Matrix The type of the matrix to multiply with.
		/// @param B The matrix to multiply with.
		/// @return The resulting matrix.
		///
		/// This function multiplies this matrix with another matrix `B`. It checks if
		/// the number of rows in `B` matches the number of columns in this matrix.
		template<typename Matrix>
		inline Matrix mul(const Matrix& B) const {

			Matrix res;
			res.resize(N, B.cols());

			if(B.rows() != K) {
				TH_MATH_ERROR("mat::transform", B.rows(), INVALID_ARGUMENT);
				algebra::mat_error(res);
				return res;
			}

			return algebra::mat_mul(res, *this, B);
		}


   		/// Overloads the `*` operator for matrix multiplication.
   		/// @tparam Matrix The type of the matrix to multiply with.
   		/// @param B The matrix to multiply with.
   		/// @return The resulting matrix.
		template<typename Matrix>
		inline auto operator*(const Matrix& B) const {

			Matrix res;
			res.resize(N, B.cols());

			if(B.rows() != K) {
				TH_MATH_ERROR("mat::transform", B.rows(), INVALID_ARGUMENT);
				algebra::mat_error(res);
				return res;
			}

			return algebra::mat_mul(res, *this, B);
		}


		/// Matrix addition.
		/// @tparam Matrix The type of the matrix to add.
		/// @param other The matrix to add to this matrix.
		/// @return Reference to this matrix after addition.
		template<typename Matrix>
		inline mat<Type, N, K>& operator+=(const Matrix& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction.
		/// @tparam Matrix The type of the matrix to subtract.
		/// @param other The matrix to subtract from this matrix.
		/// @return Reference to this matrix after subtraction.
		template<typename Matrix>
		inline mat<Type, N, K>& operator-=(const Matrix& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication.
		/// @param scalar The scalar value to multiply with.
		/// @return Reference to this matrix after multiplication.
		inline mat<Type, N, K>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division.
		/// @param scalar The scalar value to divide by.
		/// @return Reference to this matrix after division.
		///
		/// If the scalar value is close to zero, this function raises a division by zero error.
		inline mat<Type, N, K>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}

	
		/// Matrix multiplication with an assignment operator.
		/// @tparam Matrix The type of the matrix to multiply with.
		/// @param B The matrix to multiply with.
		/// @return Reference to this matrix after multiplication.
		template<typename Matrix>
		inline mat<Type, N, K>& operator*=(const Matrix& B) {
			return (*this = this->operator*(B));
		}


		/// Transposes the matrix in place.
		/// @return Reference to this matrix after transposition.
		///
		/// This function only works if the matrix is square. An assertion will trigger
		/// if the matrix is not square.
		inline mat<Type, N, K>& transpose() {
			static_assert(
				N == K, "The matrix must be square to be transposed in place.");
			return algebra::make_transposed(*this);
		}


		/// Access the element at the given row and column, by reference.
		/// This is equivalent to using (i, j).
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& get(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Get the element at the given row and column.
		/// This is equivalent to using (i, j).
		/// @param i The row index.
		/// @param j The column index.
		/// @return The element at position (i, j) as a constant reference.
		inline const Type& get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}

	
		/// Accesses the element at the given row and column with bound checking.
		/// Throws an std::out_of_range exception if the indices are out of bounds.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& at(unsigned int i, unsigned int j) {

			if (i >= rows()) {
				throw std::out_of_range(
					"The provided row index in mat::at() is out of range"
				);
			}

			if (j >= cols()) {
				throw std::out_of_range(
					"The provided column index in mat::at() is out of range"
				);
			}

			return get(i, j);
		}


		/// Accesses the element at the given row and column with bound checking.
		/// Throws an std::out_of_range exception if the indices are out of bounds.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at position (i, j).
		inline const Type& at(unsigned int i, unsigned int j) const {

			if (i >= rows()) {
				throw std::out_of_range(
					"The provided row index in mat::at() is out of range"
				);
			}

			if (j >= cols()) {
				throw std::out_of_range(
					"The provided column index in mat::at() is out of range"
				);
			}

			return get(i, j);
		}


		/// Overloads the `()` operator to access an element by reference.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& operator()(unsigned int i, unsigned int j) {
			return get(i, j);
		}


		/// Overloads the `()` operator to access an element by value.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at position (i, j).
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return get(i, j);
		}


		/// Iterator for statically allocated matrices.
		using iterator = mat_iterator<mat<Type, N, K>, Type&>;


		/// Const iterator for statically allocated matrices.
		using const_iterator = mat_iterator<const mat<Type, N, K>, const Type&>;


		/// Returns an iterator to the first element of the matrix.
		/// @return An iterator to the beginning of the matrix.
		inline auto begin() {
			return iterator(*this, 0, 0);
		}


		/// Returns an iterator to one past the last element of the matrix.
		/// @return An iterator to the end of the matrix.
		inline auto end() {
			return iterator(*this, rows(), 0);
		}

	
		/// Returns a const iterator to the first element of the matrix.
		/// @return A const iterator to the beginning of the matrix.
		inline auto begin() const {
			return const_iterator(*this, 0, 0);
		}

	
		/// Returns a const iterator to one past the last element of the matrix.
		/// @return A const iterator to the end of the matrix.
		inline auto end() const {
			return const_iterator(*this, rows(), 0);
		}

	
		/// Returns the number of rows in the matrix.
		/// @return The number of rows.
		TH_CONSTEXPR inline unsigned int rows() const {
			return N;
		}


		/// Returns the number of columns in the matrix.
		/// @return The number of columns.
		TH_CONSTEXPR inline unsigned int cols() const {
			return K;
		}


		/// Returns the total number of elements in the matrix.
		/// @return The total number of elements (rows * columns).
		inline unsigned int size() const {
			return N * K;
		}


		/// Checks whether this matrix is equal to another matrix element-wise.
		/// @tparam Matrix The type of the other matrix.
		/// @param other The matrix to compare with.
		/// @return `true` if all elements are equal, `false` otherwise.
		template<typename Matrix>
		inline bool operator==(const Matrix& other) const {
			return algebra::mat_equals(*this, other);
		}


   		/// Checks whether this matrix is not equal to another matrix element-wise.
   		/// @tparam Matrix The type of the other matrix.
   		/// @param other The matrix to compare with.
   		/// @return `true` if any elements are unequal, `false` otherwise.
		template<typename Matrix>
		inline bool operator!=(const Matrix& other) const {
			return !algebra::mat_equals(*this, other);
		}


		/// Inverts the matrix in place.
		/// @return Reference to the inverted matrix.
		///
		/// This function is only valid for square matrices.
		inline mat<Type, N, K>& invert() {
			static_assert(N == K, "The matrix must be square to be invertible.");
			return algebra::invert(*this);
		}


		/// Compatibility function to allow for allocation
		/// or resizing of dynamic matrices. Since statically
		/// allocated matrices cannot change size, this function
		/// only checks whether the target size is the same
		/// as the matrix's.
		inline mat<Type, N, K> resize(unsigned int n, unsigned int k) const {

			if(rows() != n) {
				TH_MATH_ERROR("mat::resize", n, INVALID_ARGUMENT);
			} else if(cols() != k) {
				TH_MATH_ERROR("mat::resize", k, INVALID_ARGUMENT);
			}

			return *this;
		}


#ifndef THEORETICA_NO_PRINT

			/// Converts the matrix to a string representation.
			/// @param separator Separator between elements (default is ", ").
			/// @param parenthesis Whether to enclose each row in parentheses (default is true).
			/// @return The matrix as a formatted string.
			inline std::string to_string(
				std::string separator = ", ", bool parenthesis = true) const {

				std::stringstream res;

				for (unsigned int i = 0; i < rows(); ++i) {
						
					if(parenthesis)
						res << "(";

					for (unsigned int j = 0; j < cols(); ++j) {
						
						if(j)
							res << separator;

						if(abs(get(i, j)) < MACH_EPSILON)
							res << "0";
						else
							res << get(i, j);
					}
		
					if(parenthesis)
						res << ")" << std::endl;
				}

				return res.str();
			}


			/// Convert the matrix to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Outputs the matrix to an output stream in string format.
			/// @param out The output stream.
			/// @param obj The matrix to output.
			/// @return The modified output stream.
			inline friend std::ostream& operator<<(
				std::ostream& out, const mat<Type, N, K>& obj) {
				return out << obj.to_string();
			}

#endif

	};



	/// @class mat<Type,0,0>
	/// A generic matrix with a variable number of rows and columns.
	///
	/// @param Type The type of the elements
	///
	template<typename Type>
	class mat<Type, 0, 0> {
		public:

		/// Dynamically allocated array of the elements
		std::vector<Type> data;

		/// Number of rows
		size_t row_sz {0};

		/// Number of columns
		size_t col_sz {0};


		/// Default constructor
		mat() : row_sz(0), col_sz(0) {}


   		/// Constructor that initializes a matrix with the specified number of rows and columns.
   		/// @param n Number of rows.
   		/// @param k Number of columns.
		mat(unsigned int n, unsigned int k) {
			resize(n, k);
			algebra::mat_zeroes(*this);
		}

	
		/// Copy constructor for creating a matrix from another matrix.
		/// @tparam Matrix A compatible matrix type.
		/// @param m The matrix to copy from.
		template <
			typename Matrix, enable_matrix<Matrix>
		>
		mat(const Matrix& m) {
			resize(m.rows(), m.cols());
			algebra::mat_copy(*this, m);
		}

	
		/// Constructor that initializes a matrix from a list of rows.
		/// @tparam T The type of the initializer elements (default is `Type`).
		/// @param rows Initializer list representing the matrix rows.
		template<typename T = Type>
		inline mat(const std::initializer_list<std::initializer_list<T>>& rows) {

			unsigned int i = 0;
			unsigned int j = 0;

			for (const auto& row : rows) {

				for (const auto& x : row) {
				
					if (!i && !j) {
						this->resize(rows.size(), row.size());
					}
					else if (row.size() != col_sz) {
						TH_MATH_ERROR("mat::mat", row.size(), INVALID_ARGUMENT);
						algebra::mat_error(*this);
						return;
					}

					get(i, j) = x;
					j = (j + 1) % col_sz;
				}

				i++;
			}
		}


		/// Copy assignment operator for copying from another matrix.
		/// @tparam Matrix A compatible matrix type.
   		/// @param other The matrix to copy from.
		template<typename Matrix>
		inline mat<Type>& operator=(const Matrix& other) {
			resize(other.rows(), other.cols());
			return algebra::mat_copy(*this, other);
		}


		/// Constructor that initializes a diagonal matrix with equal entries on the diagonal.
		/// @param diagonal The value for the diagonal entries.
		/// @param n Number of rows.
		/// @param k Number of columns.
		mat(Type diagonal, unsigned int n, unsigned int k) {
			
			resize(n, k);
			algebra::mat_zeroes(*this);
			const unsigned int m = min(n, k);

			for (unsigned int i = 0; i < m; ++i)
				get(i, i) = diagonal;
		}


		/// Destructor that resets the matrix size.
		~mat() {
			row_sz = 0;
			col_sz = 0;
		}


		/// Sets all elements in the matrix to zero.
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Adds two matrices element-wise.
		/// @tparam Matrix A compatible matrix type.
		/// @param other The matrix to add.
		/// @return A new matrix containing the sum of the two matrices.
		template<typename Matrix>
		inline mat<Type> operator+(const Matrix& other) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_sum(res, *this, other);
		}


	    /// Subtracts another matrix element-wise.
		/// @tparam Matrix A compatible matrix type.
		/// @param other The matrix to subtract.
		/// @return A new matrix containing the difference of the two matrices.
		template<typename Matrix>
		inline mat<Type> operator-(const Matrix& other) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_diff(res, *this, other);
		}


		/// Multiplies the matrix by a scalar.
		/// @param scalar The scalar to multiply with.
		/// @return A new matrix with each element multiplied by the scalar.
		inline mat<Type> operator*(Type scalar) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_scalmul(res, scalar, *this);
		}


		/// Friend operator to enable equations of the form
		/// (T) * (mat)
		/// @param a The scalar value.
		/// @param B The matrix to multiply.		
		inline friend mat<Type> operator*(Type a, const mat<Type>& B) {
			return B * a;
		}


		/// Friend operator to enable equations of the form
		/// (vec) * (mat) (Enables vector-matrix multiplication.)
		/// @tparam VecType The vector element type.
   		/// @tparam M The vector size.
		/// @param a The vector.
		/// @param B The matrix.
		template<typename VecType, unsigned int M>
		inline friend vec<VecType, 0> operator*(
			const vec<VecType, M>& a, const mat<Type, 0, 0>& B) {
			return B * a;
		}


		/// Divides each element in the matrix by a scalar.
		/// @param scalar The scalar to divide with.
		/// @return A new matrix with each element divided by the scalar.
		inline mat<Type> operator/(Type scalar) const {

			mat<Type> res;
			res.resize(rows(), cols());

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transforms a vector by multiplying it with the matrix.
		/// @tparam Vector The vector type.
		/// @param v The vector to transform.
		/// @return The transformed vector.
		template<typename Vector>
		inline Vector transform(const Vector& v) const {

			if(v.size() != rows()) {
				TH_MATH_ERROR("mat::transform", v.size(), INVALID_ARGUMENT);
				Vector res;
				res.resize(rows());
				return algebra::vec_error(res);
			}

			return algebra::transform(*this, v);
		}


		/// Transforms a vector by multiplying it with the matrix.
   		/// @tparam N The number of elements in the result vector.
		/// @tparam K The number of elements in the input vector.
		/// @param v The vector to transform.
		template<unsigned int N = 0, unsigned int K = 0>
		inline vec<Type, N> transform(const vec<Type, K>& v) const {
			return algebra::transform(*this, v);
		}


		/// Transforms a vector by multiplying it with the matrix.
		/// @tparam N The number of elements in the result vector.
		/// @tparam K The number of elements in the input vector.
		/// @param v The vector to transform.
		template<unsigned int N = 0, unsigned int K = 0>
		inline vec<Type, N> operator*(const vec<Type, K>& v) const {
			return transform(v);
		}


		/// Transform a vector by the matrix
		// inline vec<Type> operator*(const vec<Type>& v) const {
		// 	return transform(v);
		// }


		/// Multiplies the matrix by another matrix of the same type.
		/// @param B The matrix to multiply with.
		/// @return A new matrix containing the result of the multiplication.
		///
		/// This function performs matrix multiplication between the current matrix
		/// and another matrix `B` of the same type. The resulting matrix has the
		/// same number of rows as the current matrix and the same number of columns
		/// as `B`.
		///
		/// If the number of rows in `B` does not match the number of columns in the
		/// current matrix, an error is raised and an error matrix is returned.
		inline mat<Type> mul(const mat<Type>& B) const {

			mat<Type> res;
			res.resize(rows(), B.cols());

			if(B.rows() != cols()) {
				TH_MATH_ERROR("mat::mul", B.rows(), INVALID_ARGUMENT);
				algebra::mat_error(res);
				return res;
			}

			return algebra::mat_mul(res, *this, B);
		}


		/// Multiplies the matrix by another matrix of any compatible type.
		/// @tparam Matrix The type of the matrix to multiply with.
		/// @param B The matrix to multiply with.
		/// @return A new matrix containing the result of the multiplication.
		///
		/// This function performs matrix multiplication between the current matrix
		/// and another matrix `B`. The resulting matrix has the same number of rows
		/// as the current matrix and the same number of columns as `B`.
		///
		/// If the number of rows in `B` does not match the number of columns in the
		/// current matrix, an error is raised and an error matrix is returned.
		template<typename Matrix>
		inline Matrix mul(const Matrix& B) const {

			Matrix res;
			res.resize(rows(), B.cols());

			if(B.rows() != cols()) {
				TH_MATH_ERROR("mat::mul", B.rows(), INVALID_ARGUMENT);
				algebra::mat_error(res);
				return res;
			}

			return algebra::mat_mul(res, *this, B);
		}


		/// Matrix multiplication with any matrix type.
		/// @param B The matrix to multiply with.
		/// @return The result of multiplying the current matrix by `B`.
		///
		/// This operator overload allows the current matrix to be multiplied by
		/// another matrix `B`. The result is obtained using the `mul` function.
		template<typename Matrix>
		inline auto operator*(const Matrix& B) const {
			return mul(B);
		}


		/// Matrix addition with another matrix.
		/// @param other The matrix to add to the current matrix.
		/// @return A reference to the current matrix after addition.
		///
		/// This operator overload adds the elements of another matrix `other`
		/// to the corresponding elements of the current matrix, modifying it
		/// in place.
		template<typename Matrix>
		inline mat<Type>& operator+=(const Matrix& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction with another matrix.
		/// @param other The matrix to subtract from the current matrix.
		/// @return A reference to the current matrix after subtraction.
		///
		/// This operator overload subtracts the elements of another matrix `other`
		/// from the corresponding elements of the current matrix, modifying it
		/// in place
		template<typename Matrix>
		inline mat<Type>& operator-=(const Matrix& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication of the matrix.
		/// @param scalar The scalar value to multiply each element of the matrix by.
		/// @return A reference to the current matrix after scaling.
		///
		/// This operator overload multiplies each element of the matrix by a scalar
		/// value `scalar`, modifying the current matrix in place.
		inline mat<Type>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division of the matrix.
		/// @param scalar The scalar value to divide each element of the matrix by.
		/// @return A reference to the current matrix after division.
		///
		/// Divides each element of the matrix by `scalar`. If `scalar` is close to zero
		/// (below the defined MACH_EPSILON), an error is thrown to prevent division by zero.
		inline mat<Type>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}


		/// Matrix multiplication with any matrix type.
		/// @param B The matrix to multiply with.
		/// @return A reference to the updated matrix after multiplication.
		///
		/// This operator overload multiplies the current matrix with matrix `B`,
		/// updating its values. The result is obtained by performing matrix
		/// multiplication and storing the outcome back into the current matrix.
		template<typename Matrix>
		inline mat<Type>& operator*=(const Matrix& B) {
			return (*this = this->operator*(B));
		}


		/// Transpose the current matrix in place.
		/// @return A reference to the transposed matrix (current instance).
		///
		/// Modifies the matrix in place by transposing its elements. This operation
		/// is only valid for square matrices. For non-square matrices, use `transposed()`
		/// to obtain a new transposed matrix.
		inline mat<Type>& transpose() {
			return algebra::make_transposed(*this);
		}


		/// Access the element at the given row and column, by reference.
		/// This is equivalent to using (i, j).
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& get(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[j + i * row_sz];
#else
			return data[i + j * col_sz];
#endif
		}


		/// Get the element at the given row and column.
		/// This is equivalent to using (i, j).
		/// @param i The row index.
		/// @param j The column index.
		/// @return The element at position (i, j) as a constant reference.
		inline const Type& get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[j + i * row_sz];
#else
			return data[i + j * col_sz];
#endif
		}


		/// Access a modifiable element at a specific row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at the specified position.
		///
		/// Provides modifiable access to the matrix element at row `i` and column `j`.
		inline Type& at(unsigned int i, unsigned int j) {

			if (i >= rows()) {
				throw std::out_of_range(
					"The provided row index in mat::at() is out of range"
				);
			}

			if (j >= cols()) {
				throw std::out_of_range(
					"The provided column index in mat::at() is out of range"
				);
			}

			return get(i, j);
		}


		/// Access a constant element at a specific row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at the specified position.
		///
		/// Provides constant access to the matrix element at row `i` and column `j`.
		inline const Type& at(unsigned int i, unsigned int j) const {

			if (i >= rows()) {
				throw std::out_of_range(
					"The provided row index in mat::at() is out of range"
				);
			}

			if (j >= cols()) {
				throw std::out_of_range(
					"The provided column index in mat::at() is out of range"
				);
			}

			return get(i, j);
		}


		/// Access a modifiable element at a specific row and column using the function call operator.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at the specified position.
		///
		/// This overload allows accessing matrix elements using the function call syntax,
		/// enabling expressions like `matrix(i, j)`.
		inline Type& operator()(unsigned int i, unsigned int j) {
			return get(i, j);
		}


		/// Access a constant element at a specific row and column using the function call operator.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at the specified position.
		///
		/// This overload allows accessing matrix elements using the function call syntax
		/// in a constant context, enabling expressions like `matrix(i, j)`.
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return get(i, j);
		}


		/// Iterator for dynamically allocated matrices.
		using iterator = mat_iterator<mat<Type, 0, 0>, Type&>;


		/// Const iterator for dynamically allocated matrices.
		using const_iterator = mat_iterator<const mat<Type, 0, 0>, const Type&>;


		/// Get an iterator to the first element of the matrix.
		/// @return An iterator pointing to the first element.
		///
		/// This function provides a non-const iterator to the beginning of the matrix,
		/// allowing modification of elements.
		inline auto begin() {
			return iterator(*this, 0, 0);
		}


		/// Get an iterator to one past the last element of the matrix.
		/// @return An iterator pointing to one past the last element.
		///
		/// This function provides a non-const iterator to the end of the matrix,
		/// representing one past the last element, which is used to indicate
		/// the end of the iteration.
		inline auto end() {
			return iterator(*this, rows(), 0);
		}


		/// Get a const iterator to the first element of the matrix.
		/// @return A const iterator pointing to the first element.
		///
		/// This function provides a const iterator to the beginning of the matrix,
		/// allowing read-only access to elements.
		inline auto begin() const {
			return const_iterator(*this, 0, 0);
		}


		/// Get a const iterator to one past the last element of the matrix.
		/// @return A const iterator pointing to one past the last element.
		///
		/// This function provides a const iterator to the end of the matrix,
		/// representing one past the last element, for read-only access.
		inline auto end() const {
			return const_iterator(*this, rows(), 0);
		}


		/// Get the number of rows in the matrix.
		/// @return The number of rows in the matrix.
		TH_CONSTEXPR inline unsigned int rows() const {
			return row_sz;
		}


		/// Get the number of columns in the matrix.
		/// @return The number of columns in the matrix.
		TH_CONSTEXPR inline unsigned int cols() const {
			return col_sz;
		}


		/// Get the total number of elements of the matrix, equal to (rows * columns).
		/// @return The total number of elements in the matrix.
		inline unsigned int size() const {
			return rows() * cols();
		}


		/// Check if two matrices are equal element by element.
		/// @param other The matrix to compare with.
		/// @return True if the matrices are equal, otherwise false.
		///
		/// Compares each element of the current matrix with `other`.
		/// Returns true if all corresponding elements are equal.
		template<typename Matrix>
		inline bool operator==(const Matrix& other) const {
			return algebra::mat_equals(*this, other);
		}


		/// Check if two matrices are unequal element by element.
		/// @param other The matrix to compare with.
		/// @return True if the matrices are unequal, otherwise false.
		///
		/// Compares each element of the current matrix with `other`.
		/// Returns true if any corresponding elements are unequal.
		template<typename Matrix>
		inline bool operator!=(const Matrix& other) const {
			return !algebra::mat_equals(*this, other);
		}


		/// Invert a generic square matrix.
		/// @return A reference to the modified matrix itself, now inverted.
		///
		/// Modifies the current matrix to become its inverse. This is only defined for square, non-singular matrices.
		inline mat<Type>& invert() {
			return algebra::invert(*this);
		}


		/// Set or change the size of the matrix
		/// @param rows The number of rows
		/// @param cols The number of columns
		inline mat<Type>& resize(unsigned int rows, unsigned int cols) {

			// Do nothing if the size is already correct
			if (row_sz == rows && col_sz == cols)
				return *this;

			std::vector<Type> new_data (rows * cols);

			if (data.size()) {

				size_t min_elements = min(rows, row_sz) * min(cols, col_sz);

				// Copy the overlapping elements
				for (unsigned int i = 0; i < min_elements; ++i)
					new_data[i] = data[i];
			}

			data = new_data;
			row_sz = rows;
			col_sz = cols;

			return *this;
		}


#ifndef THEORETICA_NO_PRINT

			/// Convert the matrix to string representation
			inline std::string to_string(
				std::string separator = ", ", bool parenthesis = true) const {

				std::stringstream res;

				for (unsigned int i = 0; i < rows(); ++i) {
						
					if(parenthesis)
						res << "(";

					for (unsigned int j = 0; j < cols(); ++j) {
						
						if(j)
							res << separator;

						if(abs(get(i, j)) < MACH_EPSILON)
							res << "0";
						else
							res << get(i, j);
					}
		
					if(parenthesis)
						res << ")" << std::endl;
				}

				return res.str();
			}


			/// Convert the matrix to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the matrix in string representation to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const mat<Type>& obj) {
				return out << obj.to_string();
			}

#endif

	};

}

#endif
