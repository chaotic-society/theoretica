
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
					
					at(i, j) = x;
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
				at(i, i) = diagonal;
		}


		/// Assignment operator to copy from another matrix.
		/// @tparam Matrix The type of the matrix to copy from.
		/// @param other The matrix to copy.
		/// @return Reference to this matrix after assignment.
		template<typename Matrix>
		inline mat<Type, N, K>& operator=(const Matrix& other) {
			return algebra::mat_copy(*this, other);
		}


		/// Sets all elements of the matrix to zero.
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}

	
		/// Returns a null matrix with all elements set to zero.
		/// @return A new matrix with zero values.
		inline static mat<Type, N, K> zeroes() {
			mat<Type, N, K> res;
			algebra::mat_zeroes(res);
			return res;
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

	
		/// Returns a transposed version of the matrix.
		/// @return A new matrix that is the transposed version of this matrix.
		inline mat<Type, K, N> transposed() const {
			return algebra::transpose<mat<Type, N, K>, mat<Type, K, N>>(*this);
		}


		/// Accesses the element at the given row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Accesses the element at the given row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at position (i, j).
		inline const Type& at(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Overloads the `()` operator to access an element.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at position (i, j).
		inline Type& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Overloads the `()` operator to access an element.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at position (i, j).
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return at(i, j);
		}


		/// Gets the element at the specified row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A copy of the element at position (i, j).
		inline Type get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
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


		/// Checks if the matrix is square.
		/// @return `true` if the matrix is square (N == K), `false` otherwise.
		inline bool is_square() const {
			return algebra::is_square(*this);
		}


		/// Checks if the matrix is diagonal.
		/// @return `true` if the matrix is diagonal, `false` otherwise.
		inline bool is_diagonal() const {
			return algebra::is_diagonal(*this);
		}


		/// Checks if the matrix is symmetric.
		/// @return `true` if the matrix is symmetric, `false` otherwise.
		inline bool is_symmetric() const {
			return algebra::is_symmetric(*this);
		}


		/// Compute the density of the matrix, counting the proportion
		/// of non-zero (bigger in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of non-zero elements of the matrix.
		inline real density(real tolerance = 1E-12) {
			return algebra::density(*this, tolerance);
		}


		/// Compute the sparsity of the matrix, counting the proportion
		/// of zero (smaller in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of zero elements of the matrix.
		inline real sparsity(real tolerance = 1E-12) {
			return algebra::sparsity(*this, tolerance);
		}


		/// Computes the trace of the matrix.
		/// @return The trace (sum of diagonal elements).
		inline Type trace() {
			return algebra::trace(*this);
		}

	
		/// Computes the product of the diagonal elements.
		/// @return The product of diagonal elements.
		inline Type diagonal_product() {
			return algebra::diagonal_product(*this);
		}


		/// Computes the determinant of the matrix.
		/// @return The determinant.
		///
		/// This function is only valid for square matrices.
		inline Type det() const {
			static_assert(N == K, "The matrix must be square to compute the determinant.");
			return algebra::det(*this);
		}


		/// Computes the inverse of the matrix.
		/// @return The inverse matrix.
		///
		/// This function is only valid for square matrices.
		inline mat<Type, N, K> inverse() const {
			static_assert(N == K, "The matrix must be square to be invertible.");
			return algebra::inverse(*this);
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


		// Transformation matrices


		/// Returns the identity matrix.
		/// @return An identity matrix of the current dimensions.
		inline static mat<Type, N, K> identity() {
			return algebra::identity<mat<Type, N, K>>();
		}


		/// Returns a diagonal matrix with the specified diagonal element.
		/// @param diag The value for the diagonal elements.
		/// @return A diagonal matrix with the given diagonal value.
		inline static mat<Type, N, K> diagonal(Type diag) {
			return mat<Type, N, K>(diag);
		}

		/// Returns a 4x4 matrix for translation by the vector {x, y, z}.
		/// @tparam Vector The type of the translation vector.
		/// @param t The translation vector.
		/// @return A 4x4 translation matrix.
		template<typename Vector = vec<real, N - 1>>
		inline static mat<Type, N, K> translation(Vector&& t) {
			return algebra::translation<mat<Type, N, K>>(t);
		}

	
		/// Returns a matrix for 2D rotation by theta radians.
		/// @param theta The angle of rotation in radians.
		/// @return A 2x2 rotation matrix.
		inline static mat<Type, N, K> rotation_2d(real theta) {
			static_assert(N >= 2 && K >= 2, "The matrix must be 2x2 or bigger");
			return algebra::rotation_2d<mat<Type, N, K>>(theta);
		}


		/// Returns a matrix for 3D rotation around the x-axis.
		/// @param theta The angle of rotation in radians.
		/// @return A 3x3 rotation matrix for the x-axis.
		inline static mat<Type, N, K> rotation_3d_xaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_xaxis<mat<Type, N, K>>(theta);
		}


		/// Returns a matrix for 3D rotation around the y-axis.
		/// @param theta The angle of rotation in radians.
		/// @return A 3x3 rotation matrix for the y-axis.
		inline static mat<Type, N, K> rotation_3d_yaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_yaxis<mat<Type, N, K>>(theta);
		}


		/// Returns a matrix for 3D rotation around the z-axis.
		/// @param theta The angle of rotation in radians.
		/// @return A 3x3 rotation matrix for the z-axis.
		inline static mat<Type, N, K> rotation_3d_zaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_zaxis<mat<Type, N, K>>(theta);
		}


	    /// Returns a matrix for 3D rotation around an arbitrary axis.
		/// @tparam Vector The type of the rotation axis vector.
		/// @param theta The angle of rotation in radians.
		/// @param axis The axis vector to rotate around.
		/// @return A 3x3 rotation matrix for the given axis.
		template<typename Vector = vec<real, 3>>
		inline static mat<Type, N, K> rotation_3d(real theta, Vector&& axis) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d<mat<Type, N, K>>(theta, axis);
		}


		/// Returns a perspective projection matrix.
		/// @param left The left boundary.
		/// @param right The right boundary.
		/// @param bottom The bottom boundary.
		/// @param top The top boundary.
		/// @param near The near boundary.
		/// @param far The far boundary.
		/// @return A 4x4 perspective projection matrix.
		inline static mat<Type, N, K> perspective(
			real left, real right, real bottom,
			real top, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective<mat<Type, N, K>>(
				left, right, bottom, top, near, far);
		}

	
		/// Returns a perspective projection matrix based on field of view.
		/// @param fov The field of view angle in radians.
		/// @param aspect The aspect ratio.
		/// @param near The near boundary.
		/// @param far The far boundary.
		/// @return A 4x4 perspective projection matrix.
		inline static mat<Type, N, K> perspective_fov(
			real fov, real aspect, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective_fov<mat<Type, N, K>>(fov, aspect, near, far);
		}
	
	
		/// Returns an orthographic projection matrix.
		/// @param left The left boundary.
		/// @param right The right boundary.
		/// @param bottom The bottom boundary.
		/// @param top The top boundary.
		/// @param near The near boundary.
		/// @param far The far boundary.
		/// @return A 4x4 orthographic projection matrix.
		inline static mat<Type, N, K> ortho(
			real left, real right, real bottom, real top, real near, real far) {
			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::ortho<mat<Type, N, K>>(left, right, bottom, top, near, far);
		}


		/// Return a 4x4 transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		/// @tparam Vector1, Vector2, Vector3 Types for the camera, target, and up vectors.
		/// @param camera The camera position.
		/// @param target The target point to look at.
		/// @param up The up direction vector.
		/// @return A 4x4 look-at transformation matrix.
		template<typename Vector1, typename Vector2, typename Vector3>
		inline static mat<Type, 4, 4> look_at(
			const Vector1& camera, const Vector2& target, const Vector3& up) {
			return algebra::look_at<mat<Type, 4, 4>>(camera, target, up);
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
		/// @param n Optional parameter for number of rows.
		/// @param k Optional parameter for number of columns.
		/// @return A symplectic matrix with given dimensions.
		inline static mat<Type, N, K> symplectic(unsigned int n = 0, unsigned int k = 0) {
			static_assert(N == K && (N % 2 == 0),
				"N must equal K and they should be a multiple of 2");
			return algebra::symplectic<mat<Type, N, K>>(n, k);
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

		// Container type for storage (alias for std::vector)
		template<typename T>
		using Container = std::vector<T>;

		/// Dynamically allocated array of the elements
		Container<Container<Type>> data;

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

					at(i, j) = x;
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
				at(i, i) = diagonal;
		}


		/// Destructor that deallocates memory and resets matrix size.
		~mat() {
			row_sz = 0;
			col_sz = 0;
		}


		/// Sets all elements in the matrix to zero.
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Static method that returns a null matrix with specified dimensions.
		/// @param rows Number of rows.
		/// @param cols Number of columns.
		/// @return A matrix with all elements set to zero.
		inline static mat<Type> zeroes(unsigned int rows, unsigned int cols) {
			mat<Type> res;
			res.resize(rows, cols);
			algebra::mat_zeroes(res);
			return res;
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
				algebra::vec_error(res);
				return res;
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


		/// Return a transposed copy of the current matrix.
		/// @return A new matrix that is the transpose of the current matrix.
		///
		/// Creates a new matrix that represents the transpose of the current matrix
		/// without modifying the original matrix.
		inline mat<Type> transposed() const {
			return algebra::transpose<mat<Type>, mat<Type>>(*this);
		}


		/// Access a modifiable element at a specific row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at the specified position.
		///
		/// Provides modifiable access to the matrix element at row `i` and column `j`.
		inline Type& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access a constant element at a specific row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at the specified position.
		///
		/// Provides constant access to the matrix element at row `i` and column `j`.
		inline const Type& at(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access a modifiable element at a specific row and column using the function call operator.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A reference to the element at the specified position.
		///
		/// This overload allows accessing matrix elements using the function call syntax,
		/// enabling expressions like `matrix(i, j)`.
		inline Type& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Access a constant element at a specific row and column using the function call operator.
		/// @param i The row index.
		/// @param j The column index.
		/// @return A constant reference to the element at the specified position.
		///
		/// This overload allows accessing matrix elements using the function call syntax
		/// in a constant context, enabling expressions like `matrix(i, j)`.
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return at(i, j);
		}


		/// Get a copy of the element at a specific row and column.
		/// @param i The row index.
		/// @param j The column index.
		/// @return The value of the element at the specified position.
		///
		/// Returns a copy of the matrix element at the specified row and column indices.
		inline Type get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
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


		/// Get the total number of elements of the matrix
		/// (rows * columns)
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


		/// Determine if the matrix is square (rows == columns).
		/// @return True if the matrix is square, otherwise false.
		///
		/// Checks if the number of rows and columns are equal.
		inline bool is_square() const {
			return algebra::is_square(*this);
		}


		/// Determine if the matrix is diagonal (all non-diagonal elements are zero).
		/// @return True if the matrix is diagonal, otherwise false.
		///
		/// Checks if all elements outside the main diagonal are zero.
		inline bool is_diagonal() const {
			return algebra::is_diagonal(*this);
		}


		/// Determine if the matrix is symmetric (matrix == transpose).
		/// @return True if the matrix is symmetric, otherwise false.
		///
		/// Checks if the matrix is equal to its transpose.
		inline bool is_symmetric() const {
			return algebra::is_symmetric(*this);
		}


		/// Compute the density of the matrix, counting the proportion
		/// of non-zero (bigger in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of non-zero elements of the matrix.
		inline real density(real tolerance = 1E-12) {
			return algebra::density(*this, tolerance);
		}


		/// Compute the sparsity of the matrix, counting the proportion
		/// of zero (smaller in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of zero elements of the matrix.
		inline real sparsity(real tolerance = 1E-12) {
			return algebra::sparsity(*this, tolerance);
		}


		/// Compute the trace (sum of elements on the diagonal) of a square matrix.
		/// @return The trace of the matrix.
		///
		/// Calculates the sum of all diagonal elements in a square matrix.
		inline Type trace() {
			return algebra::trace(*this);
		}


		/// Compute the product of the diagonal elements of a square matrix.
		/// @return The product of all diagonal elements.
		///
		/// Multiplies all elements on the main diagonal to get the product.
		inline Type diagonal_product() {
			return algebra::diagonal_product(*this);
		}


		/// Compute the determinant of the matrix.
		/// @return The determinant of the matrix.
		///
		/// This function returns the determinant, which is defined only for square matrices.
		/// The determinant is a scalar value that characterizes certain properties of the matrix
		inline Type det() const {
			return algebra::det(*this);
		}


		/// Compute the inverse of a generic square matrix.
		/// @return The inverse of the matrix.
		///
		/// Returns the inverse matrix, which, when multiplied by the original matrix, yields the identity matrix.
		/// If the matrix is singular (determinant is zero), this operation is undefined.
		inline mat<Type> inverse() const {
			return algebra::inverse(*this);
		}


		/// Invert a generic square matrix.
		/// @return A reference to the modified matrix itself, now inverted.
		///
		/// Modifies the current matrix to become its inverse. This is only defined for square, non-singular matrices.
		inline mat<Type>& invert() {
			return algebra::invert(*this);
		}


		/// Set or change the size of the matrix
		/// @param n The number of rows
		/// @param k The number of columns
		inline mat<Type>& resize(unsigned int n, unsigned int k) {

			// Do nothing if the size is already correct
			if(rows() == n && cols() == k)
				return *this;

			// Distinguish between row-first and column-first allocation
#ifdef THEORETICA_ROW_FIRST
			size_t size1 = n, size2 = k;
#else
			size_t size1 = k, size2 = n;
#endif

			// The matrix must be allocated anew
			if(!data.size()) {

				data.resize(size1);
				for (unsigned int i = 0; i < size1; ++i)
					data[i].resize(size2);
				
				row_sz = n;
				col_sz = k;
				algebra::mat_zeroes(*this);

			// The matrix must be reallocated
			} else {

				std::vector<std::vector<Type>> new_data;
				new_data.resize(size1);
				for (unsigned int i = 0; i < size1; ++i)
					new_data[i].resize(size2);

				// Copy data to new memory location

				// Bounds on the region to copy
				const size_t row_bound = min(rows(), n);
				const size_t col_bound = min(cols(), n);

				for (unsigned int i = 0; i < row_bound; ++i) {
					for (unsigned int j = 0; j < col_bound; ++j) {
#ifdef THEORETICA_ROW_FIRST
						new_data[i][j] = get(i, j);
#else
						new_data[j][i] = get(i, j);
#endif
					}
				}

				// If the new matrix size is bigger,
				// zero out all remaining entries
				for (unsigned int i = 0; i < n; ++i) {
					for (unsigned int j = col_bound; j < k; ++j) {


#ifdef THEORETICA_ROW_FIRST
						new_data[i][j] = (Type) 0;
#else
						new_data[j][i] = (Type) 0;
#endif
					}
				}


				// If the new matrix size is bigger,
				// zero out all remaining entries
				for (unsigned int i = row_bound; i < n; ++i) {
					for (unsigned int j = 0; j < k; ++j) {


#ifdef THEORETICA_ROW_FIRST
						new_data[i][j] = (Type) 0;
#else
						new_data[j][i] = (Type) 0;
#endif
					}
				}

				data.clear();
				data = new_data;

				row_sz = n;
				col_sz = k;
			}

			return *this;
		}


		// Transformation matrices


		/// Get the identity matrix.
		/// @param row_sz The number of rows.
		/// @param col_sz The number of columns.
		/// @return A matrix initialized as the identity matrix.
		///
		/// Creates a matrix with ones on the main diagonal and zeros elsewhere.
		inline static mat<Type> identity(
			unsigned int row_sz, unsigned int col_sz) {

			return algebra::identity<mat<Type>>(row_sz, col_sz);
		}


		/// Get a diagonal matrix with a specified diagonal value.
		/// @param diag The value to set on the main diagonal.
		/// @param row_sz The number of rows.
		/// @param col_sz The number of columns.
		/// @return A matrix initialized as a diagonal matrix.
		inline static mat<Type> diagonal(
			Type diag, unsigned int row_sz, unsigned int col_sz) {
			return mat<Type>(diag, row_sz, col_sz);
		}


		/// Get a 4x4 matrix which translates by {x, y, z}.
		/// @param t A vector containing translation values.
		/// @return A 4x4 translation matrix.
		template<typename Vector>
		inline static mat<Type> translation(Vector&& t) {
			return algebra::translation<mat<Type>>(t);
		}


		/// Get a matrix which rotates the 2D plane by `theta` radians.
		/// @param theta The angle in radians.
		/// @return A 2D rotation matrix.
		inline static mat<Type> rotation_2d(real theta) {
			return algebra::rotation_2d<mat<Type>>(theta);
		}


		/// Get a matrix which rotates by `theta` radians around the x-axis.
		/// @param theta The angle in radians.
		/// @return A 3D rotation matrix.
		inline static mat<Type> rotation_3d_xaxis(real theta) {
			return algebra::rotation_3d_xaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates by `theta` radians around the y-axis.
		/// @param theta The angle in radians.
		/// @return A 3D rotation matrix.
		inline static mat<Type> rotation_3d_yaxis(real theta) {
			return algebra::rotation_3d_yaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates by `theta` radians around the z-axis.
		/// @param theta The angle in radians.
		/// @return A 3D rotation matrix.
		inline static mat<Type> rotation_3d_zaxis(real theta) {
			return algebra::rotation_3d_zaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates by `theta` radians around the given axis.
		/// @param theta The angle in radians.
		/// @param axis A vector representing the axis of rotation.
		/// @return A 3D rotation matrix.
		template<typename Vector = vec<real, 3>>
		inline static mat<Type> rotation_3d(real theta, Vector&& axis) {
			return algebra::rotation_3d<mat<Type>>(theta, axis);
		}

		/// Get a perspective projection matrix.
		/// @param left Left coordinate of the view volume.
		/// @param right Right coordinate of the view volume.
		/// @param bottom Bottom coordinate of the view volume.
		/// @param top Top coordinate of the view volume.
		/// @param near Near coordinate of the view volume.
		/// @param far Far coordinate of the view volume.
		/// @return A perspective projection matrix.
		inline static mat<Type> perspective(
			real left, real right, real bottom,
			real top, real near, real far) {

			return algebra::perspective<mat<Type>>(
				left, right, bottom, top, near, far);
		}

		/// Get a perspective projection matrix using field of view.
		/// @param fov Field of view angle in radians.
		/// @param aspect Aspect ratio.
		/// @param near Near clipping plane.
		/// @param far Far clipping plane.
		/// @return A perspective projection matrix.
		inline static mat<Type> perspective_fov(
			real fov, real aspect, real near, real far) {

			return algebra::perspective_fov<mat<Type>>(fov, aspect, near, far);
		}

		/// Get an orthographic projection matrix.
		/// @param left Left coordinate of the view volume.
		/// @param right Right coordinate of the view volume.
		/// @param bottom Bottom coordinate of the view volume.
		/// @param top Top coordinate of the view volume.
		/// @param near Near coordinate of the view volume.
		/// @param far Far coordinate of the view volume.
		/// @return An orthographic projection matrix.
		inline static mat<Type> ortho(
			real left, real right, real bottom, real top, real near, real far) {
			return algebra::ortho<mat<Type>>(left, right, bottom, top, near, far);
		}


		/// Return a 4x4 transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		/// @param camera Position of the camera.
		/// @param target Target point to look at.
		/// @param up Upward direction vector.
		/// @return A "look at" view transformation matrix.
		template<typename Vector1, typename Vector2, typename Vector3>
		inline static mat<Type> look_at(
			const Vector1& camera, const Vector2& target, const Vector3& up) {
			return algebra::look_at<mat<Type>>(camera, target, up);
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
		/// @param rows Number of rows in the matrix.
		/// @param cols Number of columns in the matrix.
		/// @return A symplectic matrix.
		inline static mat<Type> symplectic(unsigned int rows, unsigned int cols) {
			return algebra::symplectic<mat<Type>>(rows, cols);
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
