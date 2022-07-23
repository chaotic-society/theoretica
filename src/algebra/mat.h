
///
/// @file mat.h Matrix class and operations
///

#ifndef THEORETICA_MATRIX_H
#define THEORETICA_MATRIX_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/error.h"
#include "../core/constants.h"
#include "../core/real_analysis.h"
#include "./vec.h"


namespace theoretica {

	/// 
	/// @class mat
	/// A generic matrix with real entries
	///
	/// @param N The number of rows
	/// @param K The number of columns
	///
	template<unsigned int N, unsigned int K>
	class mat {
		public:

		/// Number of elements of the matrix
		static const unsigned int SIZE = N * K;

		/// Number of rows of the matrix
		static const unsigned int ROW_SIZE = N;

		/// Number of columns of the matrix
		static const unsigned int COL_SIZE = K;

#if defined(THEORETICA_ROW_FIRST)

		// Row-first allocation
		real data[N][K];
#else

		// Column-first allocation
		real data[K][N];
#endif

		/// Initialize a matrix to all zeroes
		inline mat() {
			make_null();
		}

		/// Initialize a matrix from another one
		inline mat(const mat<N, K>& other) {
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) = other.iget(i, j);
				}
			}
		}

		/// Initialize from an initializer list
		// inline mat(std::array<vec<K>, N> rows) {

		// 	if(rows.size() != N) {
		// 		TH_MATH_ERROR("mat::mat(std::array<vec<K>, N>)", l.size(),
		// 			INVALID_ARGUMENT);
		// 		*this = mat<N, K>(nan());
		// 		return;
		// 	}

		// 	for (unsigned int i = 0; i < N; ++i) {
		// 		for (unsigned int j = 0; j < K; ++j) {
		// 			iat(i, j) = rows[i].get(j);
		// 		}
		// 	}
		// }

		/// Copy constructor
		inline mat<N, K>& operator=(const mat<N, K>& other) {
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) = other.iget(i, j);
				}
			}
			return *this;
		}

		/// Construct a diagonal matrix
		inline mat(real diagonal) {
			make_null();
			unsigned int diag_n = min(N, K);
			for (unsigned int i = 0; i < diag_n; ++i) {
				iat(i, i) = diagonal;
			}
		}

		/// Default destructor
		inline ~mat() = default;

		/// Get the l-th column of the matrix as a vector
		inline vec<K> get_column(unsigned int l) const {
			
			vec<K> column;
			for (unsigned int i = 0; i < K; ++i)
				column.at(i) = iget(i, l);

			return column;
		}

		/// Get the l-th row of the matrix as a vector
		inline vec<N> get_row(unsigned int l) const {
			vec<N> row;
			for (unsigned int i = 0; i < N; ++i) {
				row.at(i) = iget(l, i);
			}
			return row;
		}

		/// Get the l-th row of the matrix as a vector
		inline vec<K> operator[](unsigned int l) const {
			return get_row(l);
		}

		/// Se the l-th column of the matrix from a vector
		inline void set_column(unsigned int l, const vec<N>& column) {
			for (unsigned int i = 0; i < K; ++i) {
				iat(i, l) = column.data[i];
			}
		}

		/// Se the l-th row of the matrix from a vector
		inline void set_row(unsigned int l, const vec<K>& row) {
			for (unsigned int i = 0; i < N; ++i) {
				iat(l, i) = row.data[i];
			}
		}

		/// Set all elements to zero
		inline void make_null() {
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) = 0;
				}
			}
		}

		/// Matrix addition
		inline mat<N, K> operator+(const mat<N, K>& other) const {
			mat<N, K> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iget(i, j) + other.iget(i, j);
				}
			}
			return res;
		}

		/// Matrix subtraction
		inline mat<N, K> operator-(const mat<N, K>& other) const {
			mat<N, K> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iget(i, j) - other.iget(i, j);
				}
			}
			return res;
		}

		/// Scalar multiplication
		inline mat<N, K> operator*(real scalar) const {
			mat<N, K> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iget(i, j) * scalar;
				}
			}
			return res;
		}

		/// Friend operator to enable equations of the form
		/// (T) * (mat)
		inline friend mat<N, K> operator*(real a, const mat<N, K>& B) {
			return B * a;
		}

		/// Scalar division
		inline mat<N, K> operator/(real scalar) const {

			if(scalar == 0) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return mat<N, K>(nan());
			}

			mat<N, K> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iget(i, j) / scalar;
				}
			}
			return res;
		}

		/// Transform a vector v by the matrix
		inline vec<N> transform(const vec<K>& v) const {
			vec<N> res = vec<N>();
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res[i] += iget(i, j) * v.get(j);
				}
			}
			return res;
		}

		/// Transform a vector v by the matrix
		inline vec<N> operator*(const vec<K>& v) const {
			return transform(v);
		}

		/// Matrix multiplication
		template<unsigned int M>
		inline mat<N, M> transform(const mat<K, M>& B) const {

			mat<N, M> res = mat<N, M>();
			
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < M; ++j) {
					for (unsigned int k = 0; k < K; ++k) {
						res.iat(i, j) += iget(i, k) * B.iget(k, j);
					}
				}
			}

			return res;
		}

		/// Matrix multiplication
		template<unsigned int M>
		inline mat<N, M> operator*(const mat<K, M>& B) const {
			return transform(B);
		}


		/// Matrix addition
		inline mat<N, K>& operator+=(const mat<N, K>& other) {

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) += other.iat(i, j);
				}
			}

			return *this;
		}

		/// Matrix subtraction
		inline mat<N, K>& operator-=(const mat<N, K>& other) {

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) -= other.iat(i, j);
				}
			}

			return *this;
		}

		/// Scalar multiplication
		inline mat<N, K>& operator*=(real scalar) {

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) *= scalar;
				}
			}

			return *this;
		}

		/// Scalar division
		inline mat<N, K>& operator/=(real scalar) {

			if(scalar == 0) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return mat<N, K>(nan());
			}

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) /= scalar;
				}
			}

			return *this;
		}

		/// Matrix multiplication
		inline mat<N, K>& operator*=(const mat<N, K>& B) {
			return (*this = this->operator*(B));
		}


		/// Matrix transposition
		inline void transpose() {

			if(!is_square()) {
				TH_MATH_ERROR("mat::transpose", K, IMPOSSIBLE_OPERATION);
				// Set all elements to nan ?
				return;
			}

			mat<K, N> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iat(j, i);
				}
			}

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					iat(i, j) = res.iat(i, j);
				}
			}
		}

		/// Return the transposed matrix
		inline mat<K, N> transposed() const {

			if(!is_square()) {
				TH_MATH_ERROR("mat::transposed", K, IMPOSSIBLE_OPERATION);
				return diagonal(nan());
			}

			mat<K, N> res;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					res.iat(i, j) = iget(j, i);
				}
			}
			return res;
		}


		/// Compute the dot product between v1 and v2
		/// using this matrix as the product matrix
		inline real dot(const vec<K>& v1, const vec<K>& v2) const {

			vec<N> o = transform(v2);
			real result = 0;

			for (unsigned int i = 0; i < N; ++i)
				result += v1.data[i] * o.data[i];

			return result;
		}


		/// Independent at() function.
		/// Access element at i and j index,
		/// where i is always the index on rows
		/// and j is always the index on columns.
		///
		/// This function is used inside of the library
		/// to access matrix elements independently from
		/// the specific setup of storage or notation
		/// (regardless of THEORETICA_MATRIX_LEXIC and THEORETICA_ROW_FIRST)
		inline real& iat(unsigned int i, unsigned int j) {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Access element at i and j index.
		///
		/// By default, i is the index on rows and
		/// j is the index on columns.
		/// If THEORETICA_MATRIX_LEXIC is defined,
		/// the indices will instead refer to columns
		/// and rows respectively.
		inline real& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_MATRIX_LEXIC
			return iat(j, i);
#else
			return iat(i, j);
#endif
		}


		/// Access element at indices i and j
		inline real& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Getters and setters
		inline real iget(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Getters and setters
		inline real get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_MATRIX_LEXIC
			return iget(j, i);
#else
			return iget(i, j);
#endif
		}


		/// Set the element at indices i and j
		inline void set(unsigned int i, unsigned int j, real a) {
			at(i, j) = a;
		}


		/// Check whether two matrices are equal element by element
		inline bool operator==(const mat<N, K>& other) const {

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					if(iat(i, j) != other.iat(i, j))
						return false;
				}
			}

			return true;
		}

		// Matrix types

		/// Return whether the matrix is square
		inline bool is_square() const {
			return N == K;
		}

		/// Return whether the matrix is diagonal
		inline bool is_diagonal() const {

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					if(i != j && iget(i, j) != 0)
						return false;
				}
			}
			return true;
		}

		/// Return whether the matrix is symmetric
		inline bool is_symmetric() const {

			if(!is_square())
				return false;

			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < K; ++j) {
					if(i != j && iget(i, j) != iget(j, i))
						return false;
				}
			}
			return true;
		}


		/// Compute the trace (sum of elements on the diagonal) of a square matrix
		inline real trace() {

			if(!is_square()) {
				TH_MATH_ERROR("mat::trace", K, IMPOSSIBLE_OPERATION);
				return nan();
			}

			real res = 0;

			for (unsigned int i = 0; i < N; ++i)
				res += iget(i, i);

			return res;
		}


		/// Compute the product of the diagonal elements of a square matrix
		inline real diagonal_product() {

			if(!is_square()) {
				TH_MATH_ERROR("mat::diagonal_product", K, IMPOSSIBLE_OPERATION);
				return nan();
			}

			real res = iat(0, 0);
			for (unsigned int i = 1; i < N; ++i)
				res *= iget(i, i);

			return res;
		}


		/// Return the determinant if the matrix is 2x2.
		/// @note No error checking is performed on the matrix size
		inline real det_2x2() const {
			return iget(0, 0) * iget(1, 1) - iget(1, 0) * iget(0, 1);
		}


		/// Return the determinant if the matrix is 3x3.
		/// @note No error checking is performed on the matrix size
		inline real det_3x3() const {
			return	iget(0, 0) * (iget(1, 1) * iget(2, 2) - iget(2, 1) * iget(1, 2)) +
					iget(0, 1) * (iget(1, 0) * iget(2, 2) - iget(2, 1) * iget(1, 2)) +
					iget(0, 2) * (iget(1, 0) * iget(2, 1) - iget(2, 0) * iget(1, 1));
		}


		/// Compute the determinant of the matrix using Gauss-Jordan lower triangularization
		inline real det_gj() const {

			if(!is_square()) {
				TH_MATH_ERROR("mat::det_gj", N, IMPOSSIBLE_OPERATION);
				return nan();
			}

			mat<N, K> A = mat<N, K>(*this);

			// Iterate on all columns
			for (unsigned int i = 0; i < N; ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < N; ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero.
						// The determinant does not change
						// when adding a row to another one
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < N; ++k) {
								A.iat(i, k) += A.iat(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						return 0;
					}
				}

				real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = i + 1; j < N; ++j) {

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.iat(j, i) * inv_pivot;

					// The coefficient does not change
					// when adding a linear combination
					// of a row to another
					for (unsigned int k = 0; k < N; ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
					}
				}
			}

			// The determinant of a (lower) triangular matrix
			// is the product of the elements on its diagonal
			return A.diagonal_product();
		}


		/// Compute the determinant of the matrix
		inline real det() const {

			if(!is_square()) {
				TH_MATH_ERROR("mat::det", K, IMPOSSIBLE_OPERATION);
				return nan();
			}
			
			if(N == 2)
				return det_2x2();

			if(N == 3)
				return det_3x3();

			return det_gj();
		}


		// Matrix inversion

		/// Compute the inverse of a 2x2 matrix.
		/// @note No error checking is performed on the matrix size
		inline mat<N, K> inverse_2x2() const {

			mat<N, K> B;

			// Exchange elements on the diagonal
			B.iat(0, 0) = iget(1, 1);
			B.iat(1, 1) = iget(0, 0);

			// Change sign of the other elements
			B.iat(1, 0) = -iget(1, 0);
			B.iat(0, 1) = -iget(0, 1);

			return B / B.det();
		}


		/// Compute the inverse of a generic square matrix
		inline mat<N, K> inverse() const {

			if(!is_square()) {
				TH_MATH_ERROR("mat::inverse", N, IMPOSSIBLE_OPERATION);
				return mat<N, K>(nan());
			}

			if(N == 2)
				return inverse_2x2();

			// Initialize the needed matrices
			// (A|B) is the augmented matrix
			mat<N, K> A = mat<N, K>(*this);
			mat<N, K> B = mat<N, K>::identity();

			// Gauss-Jordan elimination

			// Iterate on all columns
			for (unsigned int i = 0; i < N; ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < N; ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < N; ++k) {
								A.iat(i, k) += A.iat(j, k);
								B.iat(i, k) += B.iat(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						TH_MATH_ERROR("mat::inverse", flag, IMPOSSIBLE_OPERATION);
						return mat<N, K>(nan());
					}
				}

				real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < N; ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.iat(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < N; ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
						B.iat(j, k) -= coeff * B.iat(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < N; ++j) {
					A.iat(i, j) *= inv_pivot;
					B.iat(i, j) *= inv_pivot;
				}
				
			}

			return B;
		}


		/// Invert a generic square matrix
		inline mat<N, K>& invert() {

			if(!is_square()) {
				TH_MATH_ERROR("mat::invert", N, IMPOSSIBLE_OPERATION);
				return *this;
			}

			if(N == 2) {
				*this = inverse_2x2();
			}

			// Initialize the needed matrices
			// (A|B) is the augmented matrix
			mat<N, K> A = mat<N, K>(*this);
			mat<N, K> B = mat<N, K>::identity();

			// Gauss-Jordan elimination

			// Iterate on all columns
			for (unsigned int i = 0; i < N; ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < N; ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < N; ++k) {
								A.iat(i, k) += A.iat(j, k);
								B.iat(i, k) += B.iat(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						TH_MATH_ERROR("mat::invert", flag, IMPOSSIBLE_OPERATION);
						return *this;
					}
				}

				real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < N; ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.iat(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < N; ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
						B.iat(j, k) -= coeff * B.iat(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < N; ++j) {
					A.iat(i, j) *= inv_pivot;
					B.iat(i, j) *= inv_pivot;
				}
				
			}

			// Modify the matrix only when the inversion
			// has succeeded
			*this = B;
		}


		// Transformation matrices

		/// Get the identity matrix
		inline static mat<N, K> identity() {
			return mat<N, K>(1.0);
		}

		/// Get a diagonal matrix
		inline static mat<N, K> diagonal(real diag) {
			return mat<N, K>(diag);
		}

		/// Get a 4x4 matrix which translates by {x, y, z}
		inline static mat<4, 4> translation(vec3 t) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.iat(0, 3) = t[0];
			m.iat(1, 3) = t[1];
			m.iat(2, 3) = t[2];

			return m;
		}

		/// Get a 4x4 matrix which translates by {x, y, z}
		inline static mat<4, 4> translation(real x, real y, real z) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.iat(0, 3) = x;
			m.iat(1, 3) = y;
			m.iat(2, 3) = z;

			return m;
		}

		/// Get a matrix which rotates the 2D plane of <theta> radians
		inline static mat<2, 2> rotation_2x2(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<2, 2> res;
			res.iat(0, 0) = c;
			res.iat(0, 1) = -s;
			res.iat(1, 0) = s;
			res.iat(1, 1) = c;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<4, 4> rotation_x_4x4(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.iat(0, 0) = 1;
			res.iat(1, 1) = c;
			res.iat(2, 2) = c;
			res.iat(3, 3) = 1;

			res.iat(1, 2) = -s;
			res.iat(2, 1) = s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<3, 3> rotation_x_3x3(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.iat(0, 0) = 1;
			res.iat(1, 1) = c;
			res.iat(2, 2) = c;

			res.iat(1, 2) = -s;
			res.iat(2, 1) = s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<4, 4> rotation_y_4x4(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.iat(0, 0) = c;
			res.iat(1, 1) = 1;
			res.iat(2, 2) = c;
			res.iat(3, 3) = 1;

			res.iat(0, 2) = s;
			res.iat(2, 0) = -s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<3, 3> rotation_y_3x3(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.iat(0, 0) = c;
			res.iat(1, 1) = 1;
			res.iat(2, 2) = c;

			res.iat(0, 2) = s;
			res.iat(2, 0) = -s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<4, 4> rotation_z_4x4(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.iat(0, 0) = c;
			res.iat(1, 1) = c;
			res.iat(2, 2) = 1;
			res.iat(3, 3) = 1;

			res.iat(0, 1) = -s;
			res.iat(1, 0) = s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<3, 3> rotation_z_3x3(real theta) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.iat(0, 0) = c;
			res.iat(1, 1) = c;
			res.iat(2, 2) = 1;

			res.iat(0, 1) = -s;
			res.iat(1, 0) = s;

			return res;
		}

		/// Get a matrix which rotates <theta> radians around the <axis> axis
		inline static mat<4, 4> rotation_4x4(real theta, const vec3& axis) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			real Rx = axis.get(0);
			real Ry = axis.get(1);
			real Rz = axis.get(2);

			real cm1 = (1 - c);

			mat<4, 4> m;

			m.iat(0, 0) = c + square(Rx) * cm1;
			m.iat(0, 1) = Rx * Ry * cm1 - Rz * s;
			m.iat(0, 2) = Rx * Rz * cm1 - Ry * s;
			m.iat(0, 3) = 0;

			m.iat(1, 0) = Ry * Rx * cm1 + Rz * s;
			m.iat(1, 1) = c + square(Ry) * cm1;
			m.iat(1, 2) = Ry * Rz * cm1 - Rx * s;
			m.iat(1, 3) = 0;

			m.iat(2, 0) = Rz * Rx * cm1 - Ry * s;
			m.iat(2, 1) = Rz * Ry * cm1 + Rx * s;
			m.iat(2, 2) = c + square(Rz) * cm1;
			m.iat(2, 3) = 0;

			m.iat(3, 0) = 0;
			m.iat(3, 1) = 0;
			m.iat(3, 2) = 0;
			m.iat(3, 3) = 1;

			return m;
		}


		/// Get a matrix which rotates <theta> radians around the <axis> axis
		inline static mat<3, 3> rotation_3x3(real theta, const vec3& axis) {

			real s = theoretica::sin(theta);
			real c = theoretica::cos(theta);

			real Rx = axis.get(0);
			real Ry = axis.get(1);
			real Rz = axis.get(2);

			real cm1 = (1 - c);

			mat<3, 3> m;

			m.iat(0, 0) = c + square(Rx) * cm1;
			m.iat(0, 1) = Rx * Ry * cm1 - Rz * s;
			m.iat(0, 2) = Rx * Rz * cm1 - Ry * s;

			m.iat(1, 0) = Ry * Rx * cm1 + Rz * s;
			m.iat(1, 1) = c + square(Ry) * cm1;
			m.iat(1, 2) = Ry * Rz * cm1 - Rx * s;

			m.iat(2, 0) = Rz * Rx * cm1 - Ry * s;
			m.iat(2, 1) = Rz * Ry * cm1 + Rx * s;
			m.iat(2, 2) = c + square(Rz) * cm1;

			return m;
		}


		/// Get a scaling matrix by factors <x>, <y> and <z>
		inline static mat<4, 4> scaling_4x4(real x, real y, real z) {

			mat<4, 4> res;
			res.make_null();

			res.iat(0, 0) = x;
			res.iat(1, 1) = y;
			res.iat(2, 2) = z;
			res.iat(3, 3) = 1;

			return res;
		}


		/// Get a scaling matrix by factors <x>, <y> and <z>
		inline static mat<3, 3> scaling_3x3(real x, real y, real z) {

			mat<3, 3> res;
			res.make_null();

			res.iat(0, 0) = x;
			res.iat(1, 1) = y;
			res.iat(2, 2) = z;

			return res;
		}

		/// Get a scaling matrix by <v> factors
		template<unsigned int M>
		inline static mat<N, K> scaling(vec<M> v) {

			mat<N, K> res = mat<N, K>(1.0);

			for (unsigned int i = 0; i < min(min(M, N), K); ++i)
				res.iat(i, i) = v.get(i);

			return res;
		}


		inline static mat<4, 4> perspective(real left, real right, real bottom, real top, real near, real far) {

			mat<4, 4> result;
			result.make_null();

			result.iat(0, 0)  = 2 * near / (right - left);
			result.iat(2, 0)  = (right + left) / (right - left);
			result.iat(1, 1)  = 2 * near / (top - bottom);
			result.iat(2, 1)  = (top + bottom) / (top - bottom);
			result.iat(2, 2) = -(far + near) / (far - near);
			result.iat(3, 2) = -(2 * far * near) / (far - near);
			result.iat(2, 3) = -1;
			result.iat(3, 3) = 0;

			return result;
		}


		inline static mat<4, 4> perspective(real fov, real aspect, real near, real far) {

			real height = near * tan(radians(fov / 2.f));
			real width = height * aspect;

			return perspective(-width, width, -height, height, near, far);
		}


		inline static mat<4, 4> ortho(real left, real right, real bottom, real top, real near, real far) {

			mat<4, 4> result;
			result.make_null();

			result.iat(0, 0)  = 2 / (right - left);
			result.iat(3, 0)  = -(right + left) / (right - left);
			result.iat(1, 1)  = 2 / (top - bottom);
			result.iat(3, 1)  = -(top + bottom) / (top - bottom);
			result.iat(2, 2) = -2 / (far - near);
			result.iat(3, 2) = -(far + near) / (far - near);

			return result;
		}

		/// Return a transformation matrix that points the field of view towards a given point from the <camera> point
		inline static mat<4, 4> lookAt(vec3 camera, vec3 target, vec3 up) {

			// Construct an orthonormal basis
			vec3 z_axis = (target - camera).normalized();
			vec3 x_axis = cross(z_axis, up).normalized();
			vec3 y_axis = cross(x_axis, z_axis); // Already normalized

			// Negate z_axis to have a right-handed system
			z_axis = -z_axis;

			// Construct the rotation and translation matrix
			mat<4, 4> res;

			res.iat(0, 0) = x_axis[0];
			res.iat(0, 1) = x_axis[1];
			res.iat(0, 2) = x_axis[2];
			res.iat(0, 3) = -x_axis.dot(camera);

			res.iat(1, 0) = y_axis[0];
			res.iat(1, 1) = y_axis[1];
			res.iat(1, 2) = y_axis[2];
			res.iat(1, 3) = -y_axis.dot(camera);

			res.iat(2, 0) = z_axis[0];
			res.iat(2, 1) = z_axis[1];
			res.iat(2, 2) = z_axis[2];
			res.iat(2, 3) = -z_axis.dot(camera);

			res.iat(3, 0) = 0;
			res.iat(3, 1) = 0;
			res.iat(3, 2) = 0;
			res.iat(3, 3) = 1;

			return res;
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
		inline static mat<N, K> symplectic() {

			if(N != K || (N % 2 != 0)) {
				TH_MATH_ERROR("mat::symplectic", N, IMPOSSIBLE_OPERATION);
				return mat<N, K>(nan());
			}

			mat<N, K> res = mat<N, K>();

			for (unsigned int i = 0; i < N / 2; ++i) {	
				for (unsigned int j = N / 2; j < N; ++j) {
					if(i == (j - N / 2))
						res.iat(i, j) = 1;
				}
			}

			for (unsigned int i = N / 2; i < N; ++i) {
				for (unsigned int j = 0; j < N / 2; ++j) {
					if((i - N / 2) == (j))
						res.iat(i, j) = -1;
				}
			}

			return res;
		}



#ifndef THEORETICA_NO_PRINT

			/// Convert the matrix to string representation
			inline std::string to_string(std::string separator = ", ") const {

				std::stringstream res;

				for (unsigned int i = 0; i < N; ++i)
					res << get_row(i).to_string(separator) << std::endl;

				return res.str();
			}


			/// Stream the matrix in string representation to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const mat<N, K>& obj) {
				return out << obj.to_string();
			}

#endif

	};

	// Square matrices types

	/// A 2x2 matrix with real entries
	using mat2 = mat<2, 2>;

	/// A 3x3 matrix with real entries
	using mat3 = mat<3, 3>;

	/// A 4x4 matrix with real entries
	using mat4 = mat<4, 4>;


	/// Solve a linear system
	template<unsigned int N>
	inline vec<N> solve(mat<N, N> A, vec<N> b) {
		return A.inverse() * b;
	}


	/// Solve a linear system
	template<unsigned int N>
	inline mat<N, N> solve(mat<N, N> A, mat<N, N> B) {
		return A.inverse() * B;
	}

}

#endif
