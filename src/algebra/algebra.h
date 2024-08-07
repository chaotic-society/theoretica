
///
/// @file algebra.h Linear algebra routines.
/// This file implements all linear algebra routines of the library,
/// using templates and type traits.
/// The Matrix type must be a class with these methods:
/// - operator()()		Get the element at the i-th row and j-th column.
/// - rows()			Get the number of rows of the matrix (not defined for vectors)
/// - cols()			Get the number of columns of the matrix (not defined for vectors)
/// - resize()			Change or set the size of the matrix
///						(may throw for statically allocated matrices)
/// The Vector type must be a class with these methods:
/// - operator[]()		Get the element at index i by reference or const reference
/// - size()			Get the total number of elements of the vector
/// - resize()			Change or set the size of the vector
///						(may throw for statically allocated vectors)
///

#ifndef THEORETICA_ALGEBRA_H
#define THEORETICA_ALGEBRA_H

#include "../complex/complex_types.h"
#include "../core/core_traits.h"
#include "../core/error.h"


namespace theoretica {


	/// @namespace theoretica::algebra Linear algebra routines
	namespace algebra {


		// Operations involving one matrix or vector


		/// Overwrite the given matrix with the error
		/// matrix with NaN values on the diagonal and
		/// zeroes everywhere else. This function is used
		/// to signal an error.
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& mat_error(Matrix& m) {

			using Type = matrix_element_t<Matrix>;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m(i, j) = (Type) (i == j ? nan() : 0);

			return m;
		}


		/// Overwrite the given vector with the error
		/// vector with NaN values. This function is used
		/// to signal an error.
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_error(Vector& v) {

			using Type = indexable_element_t<Vector>;

			for (unsigned int i = 0; i < v.size(); ++i)
				v[i] = Type(nan());

			return v;
		}


		/// Overwrite a matrix with the identity matrix
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& make_identity(Matrix& m) {

			using Type = matrix_element_t<Matrix>;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m(i, j) = (Type) (i == j ? 1 : 0);

			return m;
		}


		/// Overwrite a matrix with all zeroes
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& mat_zeroes(Matrix& m) {

			using Type = matrix_element_t<Matrix>;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m(i, j) = (Type) 0;

			return m;
		}


		/// Overwrite a vector with all zeroes
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_zeroes(Vector& v) {

			using Type = indexable_element_t<Vector>;

			for (unsigned int i = 0; i < v.size(); ++i)
				v[i] = Type(0.0);

			return v;
		}


		/// Copy a matrix by overwriting another
		/// @param dest The matrix to overwrite
		/// @param src The matrix to copy
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& mat_copy(Matrix1& dest, const Matrix2& src) {

			dest.resize(src.rows(), src.cols());

			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest(i, j) = src(i, j);

			return dest;
		}


		/// Copy a vector by overwriting another.
		/// Equivalent to the operation dest = src
		/// @param dest The vector to overwrite
		/// @param src The vector to copy
		/// @return A reference to the overwritten matrix
		template<typename Vector1, typename Vector2>
		inline Vector1& vec_copy(Vector1& dest, const Vector2& src) {

			dest.resize(src.size());

			for (unsigned int i = 0; i < src.size(); ++i)
				dest[i] = src[i];

			return dest;
		}


		/// Swap two rows of a matrix, given the matrix and the
		/// two indices of the rows.
		///
		/// @param A The matrix to swap the rows of
		/// @param row1 The index of the first row to swap
		/// @param row2 The index of the other row to swap
		template<typename Matrix>
		inline Matrix& mat_swap_rows(Matrix& A, unsigned int row1, unsigned int row2) {

			using Type = matrix_element_t<Matrix>;

			if (row1 >= A.rows()) {
				TH_MATH_ERROR("algebra::mat_swap_rows", row1, INVALID_ARGUMENT);
				return mat_error(A);
			}

			if (row2 >= A.rows()) {
				TH_MATH_ERROR("algebra::mat_swap_rows", row2, INVALID_ARGUMENT);
				return mat_error(A);
			}

			if (row1 == row2)
				return A;

			for (unsigned int j = 0; j < A.cols(); ++j) {

				const Type tmp = A(row1, j);

				A(row1, j) = A(row2, j);
				A(row2, j) = tmp;
			}

			return A;
		}


		/// Swap two columns of a matrix, given the matrix and the
		/// two indices of the columns.
		///
		/// @param A The matrix to swap the columns of
		/// @param col1 The index of the first column to swap
		/// @param col2 The index of the other column to swap
		template<typename Matrix>
		inline Matrix& mat_swap_cols(Matrix& A, unsigned int col1, unsigned int col2) {

			using Type = matrix_element_t<Matrix>;

			if (col1 >= A.cols()) {
				TH_MATH_ERROR("algebra::mat_swap_cols", col1, INVALID_ARGUMENT);
				return mat_error(A);
			}

			if (col2 >= A.cols()) {
				TH_MATH_ERROR("algebra::mat_swap_cols", col2, INVALID_ARGUMENT);
				return mat_error(A);
			}

			if (col1 == col2)
				return A;

			for (unsigned int i = 0; i < A.cols(); ++i) {

				const Type tmp = A(i, col1);

				A(i, col1) = A(i, col2);
				A(i, col2) = tmp;
			}

			return A;
		}


		/// Returns the square of the Euclidean/Hermitian norm
		/// of the given vector
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto sqr_norm(const Vector& v) {

			using Type = indexable_element_t<Vector>;
			Type sum = (Type) 0;

			// Use conjugation for complex numbers
			if(is_complex_type<Type>())
				for (unsigned int i = 0; i < v.size(); ++i)
					sum += v[i] * conjugate(v[i]);
			else
				for (unsigned int i = 0; i < v.size(); ++i)
					sum += v[i] * v[i];

			return sum;
		}


		/// Returns the Euclidean/Hermitian norm
		/// of the given vector
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto norm(const Vector& v) {

			return indexable_element_t<Vector>(sqrt(sqr_norm(v)));
		}


		/// Returns the normalized vector
		/// @param v The vector to normalize
		/// @return The normalized vector
		template<typename Vector>
		inline Vector normalize(const Vector& v) {

			Vector r;
			r.resize(v.size());
			vec_copy(r, v);

			const auto m = norm(v);

			if(abs(m) < MACH_EPSILON) {
				TH_MATH_ERROR("algebra::normalize", m, DIV_BY_ZERO);
				vec_error(r);
				return r;
			}

			for (unsigned int i = 0; i < r.size(); ++i)
				r[i] /= m;

			return r;
		}


		/// Normalize a given vector overwriting it.
		/// @param v The vector to normalize
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& make_normalized(Vector& v) {

			const auto m = norm(v);

			if(abs(m) < MACH_EPSILON) {
				TH_MATH_ERROR("algebra::make_normalized", m, DIV_BY_ZERO);
				vec_error(v);
				return v;
			}

			for (unsigned int i = 0; i < v.size(); ++i)
				v[i] /= m;

			return v;
		}


		/// Computes the dot product between two vectors
		/// @param v1 The first vector
		/// @param v2 The second vector
		/// @return The dot product of the two vectors
		template<typename Vector1, typename Vector2>
		inline auto dot(const Vector1& v1, const Vector2& v2) {

			using Type = indexable_element_t<Vector1>;

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::dot", v1.size(), INVALID_ARGUMENT);
				return (Type) nan();
			}

			Type sum = 0;

			// Use conjugation for complex numbers
			if TH_CONSTIF (is_complex_type<Type>())
				for (unsigned int i = 0; i < v1.size(); ++i)
					sum += v1[i] * conjugate(v2[i]);
			else
				for (unsigned int i = 0; i < v1.size(); ++i)
					sum += v1[i] * v2[i];

			return sum;
		}


		/// Compute the cross product between two 3D vectors
		/// @param v1 The first 3D vector
		/// @param v2 The second 3D vector
		/// @return The cross product of the two vectors
		template<typename Vector1, typename Vector2>
		inline Vector1 cross(const Vector1& v1, const Vector2& v2) {

			Vector1 v3;
			v3.resize(3);

			if(v1.size() != 3) {
				TH_MATH_ERROR("algebra::cross", v1.size(), INVALID_ARGUMENT);
				vec_error(v3);
				return v3;
			}

			if(v1.size() != 3) {
				TH_MATH_ERROR("algebra::cross", v2.size(), INVALID_ARGUMENT);
				vec_error(v3);
				return v3;
			}

			v3[0] = v1[1] * v2[2] - v1[2] * v2[1];
			v3[1] = v1[2] * v2[0] - v1[0] * v2[2];
			v3[2] = v1[0] * v2[1] - v1[1] * v2[0];

			return v3;
		}


		/// Returns the transpose of the given matrix.
		/// Equivalent to the operation m^T
		/// @param m The matrix to transpose
		/// @return The transposed matrix
		template<typename Matrix, typename MatrixT = Matrix>
		inline MatrixT transpose(const Matrix& m) {

			MatrixT res;
			res.resize(m.cols(), m.rows());

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					res(j, i) = m(i, j);

			return res;
		}


		/// Transpose the given matrix
		/// @param m The matrix to transpose
		/// @return A reference to the transposed matrix
		template<typename Matrix>
		inline Matrix& make_transposed(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::make_transposed", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			for (unsigned int i = 0; i < m.rows(); ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					
					const auto buff = m(i, j);
					m(i, j) = m(j, i);
					m(j, i) = buff;
				}
			}

			return m;
		}

		
		/// Compute the transpose matrix and write the
		/// result to another matrix.
		/// Equivalent to the operation dest = src^T 
		/// @param dest The matrix to overwrite
		/// @param src The matrix to transpose
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& transpose(Matrix1& dest, const Matrix2& src) {

			// Check that the two matrices have the correct
			// number of rows and columns
			if(src.rows() != dest.cols()) {
				TH_MATH_ERROR("algebra::transpose", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.rows()) {
				TH_MATH_ERROR("algebra::transpose", src.cols(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			// Overwrite dest with the transpose of src
			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest(j, i) = src(i, j);

			return dest;
		}


		/// Returns the hermitian of the given matrix.
		/// Equivalent to the operation m^T
		/// @param m The matrix to compute the hermitian of
		/// @return The hermitian of the matrix
		template<typename Matrix, typename MatrixT = Matrix>
		inline MatrixT hermitian(const Matrix& m) {

			MatrixT res;
			res.resize(m.cols(), m.rows());

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					res(j, i) = conjugate(m(i, j));

			return res;
		}


		/// Compute the hermitian of a given matrix
		/// and overwrite it.
		/// Equivalent to the operation m = m^H
		/// @param m The matrix to compute the hermitian of
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& make_hermitian(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::hermitian", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			for (unsigned int i = 0; i < m.rows(); ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					
					const auto buff = m(i, j);
					m(i, j) = conjugate(m(j, i));
					m(j, i) = conjugate(buff);
				}
			}

			return m;
		}


		/// Hermitian (conjugate transpose) of a matrix.
		/// Equivalent to the operation dest = src^H.
		/// The base type of the matrix needs to have
		/// a compatible conjugate() function.
		/// @param dest The matrix to overwrite
		/// @param src The matrix to compute the hermitian of
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& hermitian(Matrix1& dest, const Matrix2& src) {

			// Check that the two matrices have the correct
			// number of rows and columns
			if(src.rows() != dest.cols()) {
				TH_MATH_ERROR("algebra::hermitian", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.rows()) {
				TH_MATH_ERROR("algebra::hermitian", src.cols(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			// Overwrite dest with the transpose of src
			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest(j, i) = conjugate(src(i, j));

			return dest;
		}


		/// Invert the given matrix.
		/// Equivalent to the operation dest = src^-1
		/// @param dest The matrix to overwrite
		/// @param src The matrix to invert
		/// @return A reference to the inverted matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& inverse(Matrix1& dest, const Matrix2& src) {

			if(src.rows() != src.cols()) {
				TH_MATH_ERROR("algebra::inverse", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(dest.rows() != src.rows()) {
				TH_MATH_ERROR("algebra::inverse", dest.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(dest.cols() != src.cols()) {
				TH_MATH_ERROR("algebra::inverse", dest.cols(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			using Type = matrix_element_t<Matrix2>;

			// Prepare extended matrix (A|B)
			Matrix1 A;
			A.resize(src.rows(), src.cols());
			dest.resize(src.rows(), src.cols());
			mat_copy(A, src);
			make_identity(dest);

			// Iterate on all columns
			for (unsigned int i = 0; i < src.rows(); ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A(i, i) == (Type) 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < src.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A(j, i) != (Type) 0) {

							for (unsigned int k = 0; k < src.rows(); ++k) {
								A(i, k) += A(j, k);
								dest(i, k) += dest(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						TH_MATH_ERROR("algebra::inverse", flag, IMPOSSIBLE_OPERATION);
						mat_error(dest);
						return dest;
					}
				}

				auto inv_pivot = ((Type) 1.0) / A(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < src.rows(); ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					const auto coeff = A(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < src.rows(); ++k) {
						A(j, k) -= coeff * A(i, k);
						dest(j, k) -= coeff * dest(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < src.rows(); ++j) {
					A(i, j) *= inv_pivot;
					dest(i, j) *= inv_pivot;
				}
				
			}

			return dest;
		}


		/// Returns the inverse of the given matrix.
		/// Equivalent to the operation m^-1
		/// @param m The matrix to invert
		/// @return The inverted matrix
		template<typename Matrix, typename MatrixInv = Matrix>
		inline MatrixInv inverse(const Matrix& m) {
			MatrixInv res;
			res.resize(m.rows(), m.cols());
			inverse(res, m);
			return res;
		}


		/// Invert the given matrix and overwrite it.
		/// Equivalent to the operation m = m^-1
		/// @param m The matrix to invert
		/// @return A reference to the inverted matrix
		template<typename Matrix>
		inline Matrix& invert(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::invert", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			// Prepare extended matrix (A|B)
			Matrix tmp;
			tmp.resize(m.rows(), m.cols());
			inverse(tmp, m);

			// Modify the matrix only when the inversion
			// has succeeded
			mat_copy(m, tmp);
			return m;
		}


		/// Compute the trace of the given matrix
		/// @param m A matrix of any type
		/// @return The trace of the matrix
		template<typename Matrix>
		inline auto trace(const Matrix& m) {

			auto sum = m(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				sum += m(i, i);

			return sum;
		}


		/// Compute the product of the elements of
		/// the main diagonal of a generic matrix
		/// @param m The input matrix
		/// @return The product of all the elements of the
		/// main diagonal of the input matrix
		template<typename Matrix>
		inline auto diagonal_product(const Matrix& m) {

			auto mul = m(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				mul *= m(i, i);

			return mul;
		}


		/// Compute the determinant of a square matrix.
		/// Gauss Jordan elimination is used to reduce the
		/// matrix to a triangular matrix.
		/// @param m The matrix to compute the determinant of
		/// @return The determinant of the matrix
		template<typename Matrix>
		inline auto det(const Matrix& m) {

			using Type = matrix_element_t<Matrix>;
			
			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::det", m.rows(), INVALID_ARGUMENT);
				return (Type) nan();
			}

			Matrix A;
			A.resize(m.rows(), m.cols());
			mat_copy(A, m);

			// Iterate on all columns
			for (unsigned int i = 0; i < A.rows(); ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < A.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero.
						// The determinant does not change
						// when adding a row to another one
						if(A(j, i) != (Type) 0) {

							for (unsigned int k = 0; k < A.rows(); ++k) {
								A(i, k) += (Type) A(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						return (Type) 0;
					}
				}

				const auto inv_pivot = ((Type) 1.0) / A(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = i + 1; j < A.rows(); ++j) {

					// Multiplication coefficient for
					// the elision of Ajk
					const auto coeff = (Type) A(j, i) * inv_pivot;

					// The coefficient does not change
					// when adding a linear combination
					// of a row to another
					for (unsigned int k = 0; k < A.rows(); ++k) {
						A(j, k) -= (Type) A(i, k) * coeff;
					}
				}
			}

			// The determinant of a (lower) triangular matrix
			// is the product of the elements on its diagonal
			return (Type) diagonal_product(A);
		}


		/// Return the determinant of a 2x2 matrix.
		/// @note No error checking is performed on the matrix size
		template<typename Matrix>
		inline real det_2x2(const Matrix& m) {
			return m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1);
		}


		/// Return the determinant if the matrix is 3x3.
		/// @note No error checking is performed on the matrix size
		template<typename Matrix>
		inline real det_3x3(const Matrix& m) {
			return	m(0, 0) * (m(1, 1) * m(2, 2) - m(2, 1) * m(1, 2)) -
					m(0, 1) * (m(1, 0) * m(2, 2) - m(2, 0) * m(1, 2)) +
					m(0, 2) * (m(1, 0) * m(2, 1) - m(2, 0) * m(1, 1));
		}


		/// Multiply a matrix by a scalar of any compatible type
		/// @param a A scalar value
		/// @param m The matrix to multiply
		/// @return A reference to the multiplied matrix
		template<typename Field, typename Matrix>
		inline Matrix& mat_scalmul(Field a, Matrix& m) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m(i, j) *= a;

			return m;
		}


		/// Multiply a matrix by a scalar of any compatible type
		/// which can be cast to the type of element of the output matrix.
		/// @param dest The matrix to overwrite with the result
		/// @param a A scalar value
		/// @param src The matrix to multiply
		/// @return A reference to the resulting matrix
		template<typename Field, typename Matrix1, typename Matrix2>
		inline Matrix1& mat_scalmul(Matrix1& dest, Field a, const Matrix2& src) {

			if(src.rows() != dest.rows()) {
				TH_MATH_ERROR("algebra::mat_scalmul", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.cols()) {
				TH_MATH_ERROR("algebra::mat_scalmul", src.cols(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest(i, j) = a * src(i, j);

			return dest;
		}


		/// Multiply a vector by a scalar of any compatible type
		/// @param a A scalar value
		/// @param v The vector to multiply
		/// @return A reference to the multiplied vector
		template<typename Field, typename Vector>
		inline Vector& vec_scalmul(Field a, Vector& v) {

			for (unsigned int i = 0; i < v.size(); ++i)
				v[i] *= a;

			return v;
		}


		/// Multiply a vector by a scalar of any compatible type
		/// which can be cast to the type of element of the output vector.
		/// @param dest The vector to overwrite with the result
		/// @param a A scalar value
		/// @param src The vector to multiply
		/// @return A reference to the resulting vector
		template<typename Field, typename Vector1, typename Vector2>
		inline Vector1& vec_scalmul(Vector1& dest, Field a, const Vector2& src) {

			if(src.size() != dest.size()) {
				TH_MATH_ERROR("algebra::vec_scalmul", src.size(), INVALID_ARGUMENT);
				vec_error(dest);
				return dest;
			}

			for (unsigned int i = 0; i < src.size(); ++i)
				dest[i] = a * src[i];

			return dest;
		}


		// Operations involving a matrix and a vector


		/// Apply a matrix transformation to a vector
		/// and store the result in the vector.
		/// Equivalent to the operation v = A * v
		/// @param A The matrix transformation
		/// @param v The vector to transform
		/// @return A reference to the overwritten vector
		template<typename Matrix, typename Vector>
		inline Vector& apply_transform(const Matrix& A, Vector& v) {

			if(v.size() != A.cols()) {
				TH_MATH_ERROR("algebra::apply_transform", v.size(), INVALID_ARGUMENT);
				vec_error(v);
				return v;
			}

			Vector res;
			res.resize(v.size());
			vec_zeroes(res);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res[i] += A(i, j) * v[j];

			vec_copy(res, v);
			return v;
		}


		/// Returns the matrix transformation of a vector.
		/// Equivalent to the operation A * v
		/// @param A The matrix transformation
		/// @param v The vector to transform
		/// @return The transformed vector
		template<typename Matrix, typename Vector>
		inline Vector transform(const Matrix& A, const Vector& v) {

			Vector res;
			res.resize(v.size());

			if(v.size() != A.cols()) {
				TH_MATH_ERROR("algebra::transform", v.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			vec_zeroes(res);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res[i] += A(i, j) * v[j];

			return res;
		}


		/// Apply a matrix transformation to a vector
		/// and store the result in the vector.
		/// Equivalent to the operation v = A * v
		/// @param res The matrix to overwrite with the result
		/// @param A The matrix transformation
		/// @param v The vector to transform
		/// @return A reference to the overwritten vector
		template<typename Matrix, typename Vector1, typename Vector2>
		inline Vector1& transform(Vector1& res, const Matrix& A, const Vector2& v) {

			if(v.size() != A.cols()) {
				TH_MATH_ERROR("algebra::transform", v.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			if(res.size() != v.size()) {
				TH_MATH_ERROR("algebra::transform", res.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			vec_zeroes(res);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res[i] += A(i, j) * v[j];

			return res;
		}


		// Operations involving multiple matrices or vectors


		/// Sum two matrices and store the result in the first matrix.
		/// Equivalent to the operation A = A + B
		/// @param A The first matrix to add and store the result
		/// @param B The second matrix to add
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix2& mat_sum(Matrix1& A, const Matrix2& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_sum", A.rows(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", A.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					A(i, j) = A(i, j) + B(i, j);

			return A;
		}


		/// Sum two matrices and store the result in another matrix
		/// Equivalent to the operation res = A + B
		/// @param A The first matrix to add
		/// @param B The second matrix to add
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline Matrix1& mat_sum(Matrix1& res, const Matrix2& A, const Matrix3& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_sum", A.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", A.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_sum", res.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", res.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res(i, j) = A(i, j) + B(i, j);

			return res;
		}


		/// Subtract two matrices and store the result in the first matrix.
		/// Equivalent to the operation A = A - B
		/// @param A The first matrix to store the result
		/// @param B The second matrix
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix2& mat_diff(Matrix1& A, const Matrix2& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_diff", A.rows(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", A.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					A(i, j) = A(i, j) - B(i, j);

			return A;
		}


		/// Subtract two matrices and store the result in another matrix
		/// Equivalent to the operation res = A - B
		/// @param A The first matrix
		/// @param B The second matrix
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline Matrix1& mat_diff(Matrix1& res, const Matrix2& A, const Matrix3& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_diff", A.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", A.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_diff", res.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", res.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res(i, j) = A(i, j) - B(i, j);

			return res;
		}


		/// Compute the linear combination of two matrices
		/// and store the result in the first matrix.
		/// Equivalent to the operation A = alpha * A + beta * B
		/// @param A The first matrix to combine and store the result
		/// @param B The second matrix to combine
		/// @return A reference to the overwritten matrix
		template<typename Field1, typename Matrix1, typename Field2, typename Matrix2>
		inline Matrix2& mat_lincomb(
			Field1 alpha, Matrix1& A, Field2 beta, const Matrix2& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.rows(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					A(i, j) = A(i, j) * alpha + B(i, j) * beta;

			return A;
		}


		/// Compute the linear combination of two matrices
		/// and store the result in the first matrix.
		/// Equivalent to the operation res = alpha * A + beta * B
		/// @param res The matrix to overwrite with the result
		/// @param alpha The first scalar parameter
		/// @param A The first matrix to combine
		/// @param beta The second scalar parameter
		/// @param B The second matrix to combine
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Field1, typename Matrix2,
			typename Field2, typename Matrix3>
		inline Matrix1& mat_lincomb(
			Matrix1& res, Field1 alpha, const Matrix2& A, Field2 beta, const Matrix3& B) {

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.rows(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_lincomb", res.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", res.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res(i, j) = A(i, j) * alpha + B(i, j) * beta;

			return A;
		}


		/// Multiply two matrices and store the result in the first matrix.
		/// Equivalent to the operation A = A * B
		/// @param A The first matrix to multiply and store the result
		/// @param B The second matrix to multiply
		/// @return A reference to the resulting matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& mat_mul(Matrix1& A, const Matrix2& B) {

			if(A.cols() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", A.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			if(B.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_mul", B.cols(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			Matrix1 res;
			res.resize(A.rows(), B.cols());
			mat_zeroes(res);
			
			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < B.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						res(i, j) += A(i, k) * B(k, j);

			A.resize(A.rows(), B.cols());
			mat_copy(A, res);
			return A;
		}


		/// Multiply two matrices and store the result in another matrix.
		/// Equivalent to the operation res = A * B
		/// @param res The matrix to overwrite with the result
		/// @param A The first matrix to multiply
		/// @param B The second matrix to multiply
		/// @return A reference to the resulting matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline Matrix1& mat_mul(Matrix1& res, const Matrix2& A, const Matrix3& B) {
			
			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", res.rows(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(res.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_mul", res.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(A.cols() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", A.cols(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			mat_zeroes(res);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < B.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						res(i, j) += A(i, k) * B(k, j);

			return res;
		}


		/// Checks whether two matrices are equal
		/// @param A The first matrix
		/// @param B The second matrix
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to 10 * MACH_EPSILON
		/// @return A boolean value
		template<typename Matrix1, typename Matrix2>
		inline bool mat_equals(
			const Matrix1& A, const Matrix2& B, real tolerance = 10 * MACH_EPSILON) {

			if(A.rows() != B.rows() || A.cols() != B.cols())
				return false;

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					if(abs(A(i, j) - B(i, j)) > tolerance)
						return false;

			return true;
		}


		/// Sum two vectors and store the result in the first vector.
		/// Equivalent to the operation v1 = v1 + v2
		/// @param v1 The first vector to add and store the result
		/// @param v2 The second vector to add
		/// @return A reference to the overwritten vector
		template<typename Vector1, typename Vector2>
		inline Vector2& vec_sum(Vector1& v1, const Vector2& v2) {

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::vec_sum", v1.size(), INVALID_ARGUMENT);
				vec_error(v1);
				return v1;
			}

			for (unsigned int i = 0; i < v1.size(); ++i)
				v1[i] = v1[i] + v2[i];

			return v1;
		}


		/// Sum two vectors and store the result in another vector
		/// Equivalent to the operation res = v1 + v2
		/// @param v1 The first vector to add
		/// @param v2 The second vector to add
		/// @return A reference to the overwritten vector
		template<typename Vector1, typename Vector2, typename Vector3>
		inline Vector1& vec_sum(Vector1& res, const Vector2& v1, const Vector3& v2) {

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::vec_sum", v1.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			if(res.size() != v1.size()) {
				TH_MATH_ERROR("algebra::vec_sum", res.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			for (unsigned int i = 0; i < v1.size(); ++i)
				res[i] = v1[i] + v2[i];

			return res;
		}


		/// Subtract two vectors and store the result in the first vector.
		/// Equivalent to the operation v1 = v1 - v2
		/// @param v1 The first vector to store the result
		/// @param v2 The second vector
		/// @return A reference to the overwritten vector
		template<typename Vector1, typename Vector2>
		inline Vector2& vec_diff(Vector1& v1, const Vector2& v2) {

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::vec_diff", v1.size(), INVALID_ARGUMENT);
				vec_error(v1);
				return v1;
			}

			for (unsigned int i = 0; i < v1.size(); ++i)
				v1[i] = v1[i] - v2[i];

			return v1;
		}


		/// Subtract two vectors and store the result in another vector
		/// Equivalent to the operation res = v1 - v2
		/// @param v1 The first vector
		/// @param v2 The second vector
		/// @return A reference to the overwritten vector
		template<typename Vector1, typename Vector2, typename Vector3>
		inline Vector1& vec_diff(Vector1& res, const Vector2& v1, const Vector3& v2) {

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::vec_diff", v1.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			if(res.size() != v1.size()) {
				TH_MATH_ERROR("algebra::vec_diff", res.size(), INVALID_ARGUMENT);
				vec_error(res);
				return res;
			}

			for (unsigned int i = 0; i < v1.size(); ++i)
				res[i] = v1[i] - v2[i];

			return res;
		}


		// Matrix properties


		/// Returns whether the matrix is square
		/// @param m The matrix to consider
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_square(const Matrix& m) {
			return (m.rows() == m.cols());
		}


		/// Returns whether the matrix is diagonal
		/// @param m The matrix to consider
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to 10 * MACH_EPSILON
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_diagonal(const Matrix& m, real tolerance = 10 * MACH_EPSILON) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					if(i != j && abs(m(i, j)) > tolerance)
						return false;

			return true;
		}

		/// Returns whether the matrix is symmetric
		/// @param m The matrix to consider
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to MATRIX_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_symmetric(const Matrix& m, real tolerance = MATRIX_ELEMENT_TOL) {

			if(!is_square(m))
				return false;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					if(abs(m(i, j) - m(j, i)) > tolerance)
						return false;

			return true;
		}


		/// Returns whether the matrix is lower triangular
		/// @param m The matrix to consider
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to MATRIX_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_lower_triangular(const Matrix& m, real tolerance = MATRIX_ELEMENT_TOL) {

			if(!is_square(m))
				return false;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = i + 1; j < m.cols(); ++j)
					if (abs(m(i, j)) > tolerance)
						return false;

			return true;
		}


		/// Returns whether the matrix is upper triangular
		/// @param m The matrix to consider
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to MATRIX_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_upper_triangular(const Matrix& m, real tolerance = MATRIX_ELEMENT_TOL) {

			if(!is_square(m))
				return false;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < i; ++j)
					if (abs(m(i, j)) > tolerance)
						return false;

			return true;
		}


		// Matrix decompositions


		/// Decompose a square matrix to two triangular matrices,
		/// L and U where L is lower and U is upper, so that \f$A = LU\f$.
		///
		/// @param A The matrix to decompose
		/// @param L The matrix to overwrite with the lower triangular one
		/// @param U The matrix to overwrite with the upper triangular one
		/// @return The LU decomposition of the matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline void decompose_lu(const Matrix1& A, Matrix2& L, Matrix3& U) {

			if (!is_square(A)) {
				TH_MATH_ERROR("algebra::decompose_lu", A.rows(), INVALID_ARGUMENT);
				mat_error(L); mat_error(U);
				return;
			}

			if (A.rows() != L.rows()) {
				TH_MATH_ERROR("algebra::decompose_lu", L.rows(), INVALID_ARGUMENT);
				mat_error(L); mat_error(U);
				return;
			}

			if (A.cols() != L.cols()) {
				TH_MATH_ERROR("algebra::decompose_lu", L.cols(), INVALID_ARGUMENT);
				mat_error(L); mat_error(U);
				return;
			}

			if (A.rows() != U.rows()) {
				TH_MATH_ERROR("algebra::decompose_lu", U.rows(), INVALID_ARGUMENT);
				mat_error(L); mat_error(U);
				return;
			}

			if (A.cols() != U.cols()) {
				TH_MATH_ERROR("algebra::decompose_lu", U.cols(), INVALID_ARGUMENT);
				mat_error(L); mat_error(U);
				return;
			}

			using Type = matrix_element_t<Matrix1>;

			// Set the diagonal of L to 1.0
			for (unsigned int i = 0; i < A.rows(); ++i)
				L(i, i) = (Type) 1.0;

			// Compute L and U
			for(unsigned int i = 0; i < A.rows(); ++i) {
				
				for(unsigned int j = 0; j < A.rows(); ++j) {

					U(i, j) = A(i, j);
				
					for(unsigned int k = 0; k < i; ++k)
						U(i, j) -= L(i, k) * U(k, j);
				}

				for(unsigned int j = i + 1; j < A.rows(); ++j) {

					L(j, i) = A(j, i);

					for(unsigned int k = 0; k < i; ++k)
						L(j, i) -= L(j, k) * U(k, i);

					L(j, i) /= U(i, i);
				}
			}
		}


		/// Decompose a square matrix to two triangular matrices,
		/// L and U where L is lower and U is upper, so that \f$A = LU\f$
		/// overwriting the matrix A with the elements of both matrices,
		/// omitting the diagonal of L (equal to all ones).
		/// Particularly useful for solving linear systems.
		///
		/// @param A The matrix to decompose and overwrite
		/// @return A reference to the overwritten matrix A
		template<typename Matrix>
		inline Matrix& decompose_lu_inplace(Matrix& A) {

			if (!is_square(A)) {
				TH_MATH_ERROR("algebra::decompose_lu_inplace", A.rows(), INVALID_ARGUMENT);
				mat_error(A);
				return A;
			}

			for(unsigned int j = 0; j < A.cols(); ++j) {

				for(unsigned int i = j + 1; i < A.rows(); ++i) {

					A(i, j) /= A(j, j);

					for(unsigned int k = j + 1; k < A.rows(); ++k)
						A(i, k) -= A(i, j) * A(j, k);
				}
			}
		
			return A;
		}


		/// Decompose a symmetric positive definite matrix into
		/// a triangular matrix so that \f$A = L L^T\f$ using
		/// Cholesky decomposition.
		///
		/// @param A The matrix to decompose
		/// @return The Cholesky decomposition of the matrix
		template<typename Matrix>
		inline Matrix decompose_cholesky(const Matrix& A) {

			Matrix L;
			L.resize(A.rows(), A.cols());
			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("algebra::decompose_cholesky", A.rows(), INVALID_ARGUMENT);
				return mat_error(L);
			}

			if (!is_symmetric(A)) {
				TH_MATH_ERROR("algebra::decompose_cholesky", false, INVALID_ARGUMENT);
				return mat_error(L);
			}

			mat_zeroes(L);

			for (unsigned int i = 0; i < A.cols(); ++i) {
				
				for (unsigned int j = 0; j <= i; ++j) {
					
					Type sum = L(i, 0) * L(j, 0);

					for (unsigned int k = 1; k < j; ++k)
						sum += L(i, k) * L(j, k);

					if (i == j) {

						const Type sqr_diag = A(j, j) - sum;

						// Additional check to ensure that the matrix is positive definite
						if (sqr_diag < MACH_EPSILON) {
							TH_MATH_ERROR("algebra::decompose_cholesky", sqr_diag, INVALID_ARGUMENT);
							return mat_error(L);
						}

						L(i, j) = sqrt(sqr_diag);

					} else {

						L(i, j) = (A(i, j) - sum) / L(j, j);
					}
				}
			}

			return L;
		}


		// Linear system solvers


		/// Solve the linear system \f$Lx = b\f$ for lower triangular \f$L\f$.
		/// @note No check is performed on the triangularity of \f$L\f$.
		///
		/// @param L The lower triangular matrix.
		/// @param b The known vector.
		template<typename Matrix, typename Vector>
		inline Vector solve_triangular_lower(const Matrix& L, const Vector& b) {

			Vector x;
			x.resize(L.cols());
			using Type = matrix_element_t<Matrix>;

			if (!is_square(L)) {
				TH_MATH_ERROR("solve_triangular_lower", false, INVALID_ARGUMENT);
				return vec_error(x);
			}

			if (b.size() != L.rows()) {
				TH_MATH_ERROR("solve_triangular_lower", b.size(), INVALID_ARGUMENT);
				return vec_error(x);
			}

			// Solve using forward substitution
			for (unsigned int i = 0; i < L.cols(); ++i) {
				
				Type sum = L(i, 0) * x[0];

				for (unsigned int j = 1; j < i; ++j)
					sum += L(i, j) * x[j];

				if (abs(L(i, i)) < MACH_EPSILON) {
					TH_MATH_ERROR("solve_triangular_lower", L(i, i), DIV_BY_ZERO);
					return vec_error(x);
				}

				x[i] = (b[i] - sum) / L(i, i);
			}

			return x;
		}


		/// Solve the linear system \f$Ux = b\f$ for upper triangular \f$U\f$.
		/// @note No check is performed on the triangularity of \f$U\f$.
		///
		/// @param U The upper triangular matrix.
		/// @param b The known vector.
		template<typename Matrix, typename Vector>
		inline Vector solve_triangular_upper(const Matrix& U, const Vector& b) {

			Vector x;
			x.resize(U.cols());
			using Type = matrix_element_t<Matrix>;

			if (!is_square(U)) {
				TH_MATH_ERROR("solve_triangular_upper", false, INVALID_ARGUMENT);
				return vec_error(x);
			}

			if (b.size() != U.rows()) {
				TH_MATH_ERROR("solve_triangular_upper", b.size(), INVALID_ARGUMENT);
				return vec_error(x);
			}

			if (U.rows() == 1) {
				x[0] = b[0] / U(0, 0);
				return x;
			}

			// Solve using backward substitution
			for (int i = U.rows() - 1; i >= 0; --i) {
				
				Type sum = U(i, i + 1) * x[i + 1];

				for (unsigned int j = i + 2; j < U.cols(); ++j)
					sum += U(i, j) * x[j];

				if (abs(U(i, i)) < MACH_EPSILON) {
					TH_MATH_ERROR("solve_triangular_upper", U(i, i), DIV_BY_ZERO);
					return vec_error(x);
				}

				x[i] = (b[i] - sum) / U(i, i);
			}

			return x;
		}


		/// Solve the linear system \f$Tx = b\f$ for triangular \f$T\f$.
		/// The correct solver is selected depending on the elements of \f$T\f$,
		/// if the property of the matrix is known a priori, calling the
		/// specific function is more efficient.
		/// @note No check is performed on the triangularity of \f$T\f$.
		///
		/// @param T The triangular matrix.
		/// @param b The known vector.
		template<typename Matrix, typename Vector>
		inline Vector solve_triangular(const Matrix& T, const Vector& b) {

			// Pick the correct solver
			if (is_lower_triangular(T))
				return solve_triangular_lower(T, b);
			else if(is_upper_triangular(T))
				return solve_triangular_upper(T, b);
			else {
				Vector err;
				err.resize(b.size());
				return vec_error(err);
			}
		}
	}
}

#endif
