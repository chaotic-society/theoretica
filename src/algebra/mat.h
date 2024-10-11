
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

			/// Construct iterator from a matrix
			mat_iterator(
				Matrix& matrix,
				size_t row = 0,
				size_t col = 0)
			: matrix(matrix), row(row), col(col) {}


			/// Dereference the iterator to get
			/// the current element by reference.
			ReturnType operator*() {
				return matrix(row, col);
			}


			/// Move to the next element in the matrix.
			mat_iterator& operator++() {

				col++;

				if(col == matrix.cols()) {
					col = 0;
					row++;
				}

				return *this;
			}


			/// Get the index of the current row
			size_t row_index() {
				return row;
			}


			/// Get the index of the current column
			size_t col_index() {
				return col;
			}


			/// Move to the previous element in the matrix.
			// mat_iterator& operator--() {
			// 	if(row == 0) {
			// 		TH_MATH_ERROR(
			// 			"iterator::operator--",
			// 			row, IMPOSSIBLE_OPERATION);
			// 		return *this;
			// 	}
			// 	if(col == 0) {
			// 		col = max_col;
			// 		row--;
			// 	}
			// 	col--;
			// 	return *this;
			// }


			// Comparison operators

			bool operator==(const mat_iterator& other) const {
				return (row == other.row) &&
					(col == other.col);
			}

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


		/// Default constructor
		mat() {
			algebra::mat_zeroes(*this);
		}


		/// Copy constructor
		template<typename Matrix>
		mat(const Matrix& m) {
			algebra::mat_copy(*this, m);
		}


		/// Construct from a list of the rows
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


		/// Copy constructor
		template<typename Matrix>
		inline mat<Type, N, K>& operator=(const Matrix& other) {
			return algebra::mat_copy(*this, other);
		}


		/// Set all elements to zero
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Get the null matrix
		inline static mat<Type, N, K> zeroes() {
			mat<Type, N, K> res;
			algebra::mat_zeroes(res);
			return res;
		}


		/// Matrix addition
		template<typename Matrix>
		inline mat<Type, N, K> operator+(const Matrix& other) const {
			mat<Type, N, K> res;
			return algebra::mat_sum(res, *this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<Type, N, K> operator-(const Matrix& other) const {
			mat<Type, N, K> res;
			return algebra::mat_diff(res, *this, other);
		}


		/// Scalar multiplication
		inline mat<Type, N, K> operator*(Type scalar) const {
			mat<Type, N, K> res;
			return algebra::mat_scalmul(res, scalar, *this);
		}


		/// Friend operator to enable equations of the form
		/// (T) * (mat)
		inline friend mat<Type, N, K> operator*(Type a, const mat<Type, N, K>& B) {
			return B * a;
		}


		/// Friend operator to enable equations of the form
		/// (vec) * (mat)
		template<typename VecType, unsigned int M>
		inline friend vec<VecType, K> operator*(
			const vec<VecType, M>& a, const mat<Type, N, K>& B) {
			return B * a;
		}


		/// Scalar division
		inline mat<Type, N, K> operator/(Type scalar) const {

			mat<Type, N, K> res;

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transform a vector v by the matrix
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


		/// Transform a vector by the matrix
		inline vec<Type, N> transform(const vec<Type, K>& v) const {
			return algebra::transform(*this, v);
		}


		/// Transform a vector by the matrix
		inline vec<Type, N> operator*(const vec<Type, K>& v) const {
			return transform(v);
		}


		/// Matrix multiplication
		template<unsigned int M>
		inline mat<Type, N, M> mul(const mat<Type, K, M>& B) const {
			mat<Type, N, M> res;
			return algebra::mat_mul(res, *this, B);
		}


		/// Matrix multiplication by any matrix type
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


		/// Matrix multiplication
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


		/// Matrix addition
		template<typename Matrix>
		inline mat<Type, N, K>& operator+=(const Matrix& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<Type, N, K>& operator-=(const Matrix& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication
		inline mat<Type, N, K>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division
		inline mat<Type, N, K>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}


		/// Matrix multiplication
		template<typename Matrix>
		inline mat<Type, N, K>& operator*=(const Matrix& B) {
			return (*this = this->operator*(B));
		}


		/// Transpose the matrix itself
		inline mat<Type, N, K>& transpose() {
			static_assert(
				N == K, "The matrix must be square to be transposed in place.");
			return algebra::make_transposed(*this);
		}


		/// Return the transposed matrix, without modifying the
		/// matrix itself.
		inline mat<Type, K, N> transposed() const {
			return algebra::transpose<mat<Type, N, K>, mat<Type, K, N>>(*this);
		}


		/// Access the element at the i-th row and j-th column
		inline Type& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access the element at the i-th row and j-th column
		inline const Type& at(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access the element at the i-th row and j-th column
		inline Type& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Get the element at the i-th row and j-th column
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return at(i, j);
		}


		/// Get the element at the i-th row and j-th column
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


		/// Get an iterator to the first element of the matrix.
		inline auto begin() {
			return iterator(*this, 0, 0);
		}


		/// Get an iterator to one plus the last element of the matrix.
		inline auto end() {
			return iterator(*this, rows(), 0);
		}


		/// Get a const iterator to the first element
		/// of the matrix.
		inline auto begin() const {
			return const_iterator(*this, 0, 0);
		}


		/// Get a const iterator to one plus the last element
		/// of the matrix.
		inline auto end() const {
			return const_iterator(*this, rows(), 0);
		}


		/// Get the number of rows of the matrix
		TH_CONSTEXPR inline unsigned int rows() const {
			return N;
		}


		/// Get the number of columns of the matrix
		TH_CONSTEXPR inline unsigned int cols() const {
			return K;
		}


		/// Get the total number of elements of the matrix
		/// (rows * columns)
		inline unsigned int size() const {
			return N * K;
		}


		/// Check whether two matrices are equal element by element
		template<typename Matrix>
		inline bool operator==(const Matrix& other) const {
			return algebra::mat_equals(*this, other);
		}


		/// Check whether two matrices are unequal element by element
		template<typename Matrix>
		inline bool operator!=(const Matrix& other) const {
			return !algebra::mat_equals(*this, other);
		}


		/// Return whether the matrix is square
		inline bool is_square() const {
			return algebra::is_square(*this);
		}


		/// Return whether the matrix is diagonal
		inline bool is_diagonal() const {
			return algebra::is_diagonal(*this);
		}


		/// Return whether the matrix is symmetric
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


		/// Compute the trace (sum of elements on the diagonal) of a square matrix
		inline Type trace() {
			return algebra::trace(*this);
		}


		/// Compute the product of the diagonal elements of a square matrix
		inline Type diagonal_product() {
			return algebra::diagonal_product(*this);
		}


		/// Compute the determinant of the matrix
		inline Type det() const {
			static_assert(N == K, "The matrix must be square to compute the determinant.");
			return algebra::det(*this);
		}


		/// Compute the inverse of a generic square matrix
		inline mat<Type, N, K> inverse() const {
			static_assert(N == K, "The matrix must be square to be invertible.");
			return algebra::inverse(*this);
		}


		/// Invert a generic square matrix
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


		/// Get the identity matrix
		inline static mat<Type, N, K> identity() {
			return algebra::identity<mat<Type, N, K>>();
		}


		/// Get a diagonal matrix
		inline static mat<Type, N, K> diagonal(Type diag) {
			return mat<Type, N, K>(diag);
		}


		/// Get a 4x4 matrix which translates by {x, y, z}
		template<typename Vector = vec<real, N - 1>>
		inline static mat<Type, N, K> translation(Vector&& t) {
			return algebra::translation<mat<Type, N, K>>(t);
		}


		/// Get a matrix which rotates the 2D plane of <theta> radians
		inline static mat<Type, N, K> rotation_2d(real theta) {
			static_assert(N >= 2 && K >= 2, "The matrix must be 2x2 or bigger");
			return algebra::rotation_2d<mat<Type, N, K>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<Type, N, K> rotation_3d_xaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_xaxis<mat<Type, N, K>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<Type, N, K> rotation_3d_yaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_yaxis<mat<Type, N, K>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<Type, N, K> rotation_3d_zaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_zaxis<mat<Type, N, K>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the <axis> axis
		template<typename Vector = vec<real, 3>>
		inline static mat<Type, N, K> rotation_3d(real theta, Vector&& axis) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d<mat<Type, N, K>>(theta, axis);
		}


		inline static mat<Type, N, K> perspective(
			real left, real right, real bottom,
			real top, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective<mat<Type, N, K>>(
				left, right, bottom, top, near, far);
		}


		inline static mat<Type, N, K> perspective_fov(
			real fov, real aspect, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective_fov<mat<Type, N, K>>(fov, aspect, near, far);
		}


		inline static mat<Type, N, K> ortho(
			real left, real right, real bottom, real top, real near, real far) {
			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::ortho<mat<Type, N, K>>(left, right, bottom, top, near, far);
		}


		/// Return a 4x4 transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		template<typename Vector1, typename Vector2, typename Vector3>
		inline static mat<Type, 4, 4> look_at(
			const Vector1& camera, const Vector2& target, const Vector3& up) {
			return algebra::look_at<mat<Type, 4, 4>>(camera, target, up);
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
		inline static mat<Type, N, K> symplectic(unsigned int n = 0, unsigned int k = 0) {
			static_assert(N == K && (N % 2 == 0),
				"N must equal K and they should be a multiple of 2");
			return algebra::symplectic<mat<Type, N, K>>(n, k);
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


		/// Construct a matrix with n rows and k columns
		mat(unsigned int n, unsigned int k) {
			resize(n, k);
			algebra::mat_zeroes(*this);
		}


		/// Copy constructor
		template <
			typename Matrix, enable_matrix<Matrix>
		>
		mat(const Matrix& m) {
			resize(m.rows(), m.cols());
			algebra::mat_copy(*this, m);
		}


		/// Construct from a list of the rows
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


		/// Copy constructor
		template<typename Matrix>
		inline mat<Type>& operator=(const Matrix& other) {
			resize(other.rows(), other.cols());
			return algebra::mat_copy(*this, other);
		}


		/// Construct a diagonal matrix with all equal entries
		/// with n rows and k columns.
		mat(Type diagonal, unsigned int n, unsigned int k) {
			
			resize(n, k);
			algebra::mat_zeroes(*this);
			const unsigned int m = min(n, k);

			for (unsigned int i = 0; i < m; ++i)
				at(i, i) = diagonal;
		}


		/// Deallocate memory
		~mat() {
			row_sz = 0;
			col_sz = 0;
		}


		/// Set all elements to zero
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Get the null matrix
		inline static mat<Type> zeroes(unsigned int rows, unsigned int cols) {
			mat<Type> res;
			res.resize(rows, cols);
			algebra::mat_zeroes(res);
			return res;
		}


		/// Matrix addition
		template<typename Matrix>
		inline mat<Type> operator+(const Matrix& other) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_sum(res, *this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<Type> operator-(const Matrix& other) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_diff(res, *this, other);
		}


		/// Scalar multiplication
		inline mat<Type> operator*(Type scalar) const {
			mat<Type> res;
			res.resize(rows(), cols());
			return algebra::mat_scalmul(res, scalar, *this);
		}


		/// Friend operator to enable equations of the form
		/// (T) * (mat)
		inline friend mat<Type> operator*(Type a, const mat<Type>& B) {
			return B * a;
		}


		/// Friend operator to enable equations of the form
		/// (vec) * (mat)
		template<typename VecType, unsigned int M>
		inline friend vec<VecType, 0> operator*(
			const vec<VecType, M>& a, const mat<Type, 0, 0>& B) {
			return B * a;
		}


		/// Scalar division
		inline mat<Type> operator/(Type scalar) const {

			mat<Type> res;
			res.resize(rows(), cols());

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transform a vector v by the matrix
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


		/// Transform a vector by the matrix
		template<unsigned int N = 0, unsigned int K = 0>
		inline vec<Type, N> transform(const vec<Type, K>& v) const {
			return algebra::transform(*this, v);
		}


		/// Transform a vector by the matrix
		template<unsigned int N = 0, unsigned int K = 0>
		inline vec<Type, N> operator*(const vec<Type, K>& v) const {
			return transform(v);
		}


		/// Transform a vector by the matrix
		// inline vec<Type> operator*(const vec<Type>& v) const {
		// 	return transform(v);
		// }


		/// Matrix multiplication
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


		/// Matrix multiplication by any matrix type
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


		/// Matrix multiplication
		template<typename Matrix>
		inline auto operator*(const Matrix& B) const {
			return mul(B);
		}


		/// Matrix addition
		template<typename Matrix>
		inline mat<Type>& operator+=(const Matrix& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<Type>& operator-=(const Matrix& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication
		inline mat<Type>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division
		inline mat<Type>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}


		/// Matrix multiplication
		template<typename Matrix>
		inline mat<Type>& operator*=(const Matrix& B) {
			return (*this = this->operator*(B));
		}


		/// Transpose the matrix itself
		inline mat<Type>& transpose() {
			return algebra::make_transposed(*this);
		}


		/// Return the transposed matrix, without modifying the
		/// matrix itself.
		inline mat<Type> transposed() const {
			return algebra::transpose<mat<Type>, mat<Type>>(*this);
		}


		/// Access the element at the i-th row and j-th column
		inline Type& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access the element at the i-th row and j-th column
		inline const Type& at(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access the element at the i-th row and j-th column
		inline Type& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Get the element at the i-th row and j-th column
		inline const Type& operator()(unsigned int i, unsigned int j) const {
			return at(i, j);
		}


		/// Get the element at the i-th row and j-th column
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
		inline auto begin() {
			return iterator(*this, 0, 0);
		}


		/// Get an iterator to one plus the last element of the matrix.
		inline auto end() {
			return iterator(*this, rows(), 0);
		}


		/// Get a const iterator to the first element of the matrix.
		inline auto begin() const {
			return const_iterator(*this, 0, 0);
		}


		/// Get a const iterator to one plus the last element of the matrix.
		inline auto end() const {
			return const_iterator(*this, rows(), 0);
		}


		/// Get the number of rows of the matrix
		TH_CONSTEXPR inline unsigned int rows() const {
			return row_sz;
		}


		/// Get the number of columns of the matrix
		TH_CONSTEXPR inline unsigned int cols() const {
			return col_sz;
		}


		/// Get the total number of elements of the matrix
		/// (rows * columns)
		inline unsigned int size() const {
			return rows() * cols();
		}


		/// Check whether two matrices are equal element by element
		template<typename Matrix>
		inline bool operator==(const Matrix& other) const {
			return algebra::mat_equals(*this, other);
		}


		/// Check whether two matrices are unequal element by element
		template<typename Matrix>
		inline bool operator!=(const Matrix& other) const {
			return !algebra::mat_equals(*this, other);
		}


		/// Return whether the matrix is square
		inline bool is_square() const {
			return algebra::is_square(*this);
		}


		/// Return whether the matrix is diagonal
		inline bool is_diagonal() const {
			return algebra::is_diagonal(*this);
		}


		/// Return whether the matrix is symmetric
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


		/// Compute the trace (sum of elements on the diagonal) of a square matrix
		inline Type trace() {
			return algebra::trace(*this);
		}


		/// Compute the product of the diagonal elements of a square matrix
		inline Type diagonal_product() {
			return algebra::diagonal_product(*this);
		}


		/// Compute the determinant of the matrix
		inline Type det() const {
			return algebra::det(*this);
		}


		/// Compute the inverse of a generic square matrix
		inline mat<Type> inverse() const {
			return algebra::inverse(*this);
		}


		/// Invert a generic square matrix
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

			// Distinguish between row-first and column-first
			// allocation
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


		/// Get the identity matrix
		inline static mat<Type> identity(
			unsigned int row_sz, unsigned int col_sz) {

			return algebra::identity<mat<Type>>(row_sz, col_sz);
		}


		/// Get a diagonal matrix
		inline static mat<Type> diagonal(
			Type diag, unsigned int row_sz, unsigned int col_sz) {
			return mat<Type>(diag, row_sz, col_sz);
		}


		/// Get a 4x4 matrix which translates by {x, y, z}
		template<typename Vector>
		inline static mat<Type> translation(Vector&& t) {
			return algebra::translation<mat<Type>>(t);
		}


		/// Get a matrix which rotates the 2D plane of <theta> radians
		inline static mat<Type> rotation_2d(real theta) {
			return algebra::rotation_2d<mat<Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<Type> rotation_3d_xaxis(real theta) {
			return algebra::rotation_3d_xaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<Type> rotation_3d_yaxis(real theta) {
			return algebra::rotation_3d_yaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<Type> rotation_3d_zaxis(real theta) {
			return algebra::rotation_3d_zaxis<mat<Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the <axis> axis
		template<typename Vector = vec<real, 3>>
		inline static mat<Type> rotation_3d(real theta, Vector&& axis) {
			return algebra::rotation_3d<mat<Type>>(theta, axis);
		}


		inline static mat<Type> perspective(
			real left, real right, real bottom,
			real top, real near, real far) {

			return algebra::perspective<mat<Type>>(
				left, right, bottom, top, near, far);
		}


		inline static mat<Type> perspective_fov(
			real fov, real aspect, real near, real far) {

			return algebra::perspective_fov<mat<Type>>(fov, aspect, near, far);
		}


		inline static mat<Type> ortho(
			real left, real right, real bottom, real top, real near, real far) {
			return algebra::ortho<mat<Type>>(left, right, bottom, top, near, far);
		}


		/// Return a 4x4 transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		template<typename Vector1, typename Vector2, typename Vector3>
		inline static mat<Type> look_at(
			const Vector1& camera, const Vector2& target, const Vector3& up) {
			return algebra::look_at<mat<Type>>(camera, target, up);
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
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
