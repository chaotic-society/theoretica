
///
/// @file algebra.h Linear algebra routines.
/// This file implements all linear algebra routines of the library,
/// working on templated data structures. The Matrix template must 
/// be a class with these methods:
/// - iget(i, j)  Get the element in the i-th row and j-th column
/// - iat(i, j)  Get a reference to the element in the i-th row and j-th column
/// - size()  Get the number of elements of the matrix or vector (N. rows * N. columns)
/// - rows()  Get the number of rows of the matrix (not defined for vectors)
/// - cols()  Get the number of columns of the matrix (not defined for vectors)
/// - resize() Change or set the size of the matrix
///

#ifndef THEORETICA_ALGEBRA_H
#define THEORETICA_ALGEBRA_H


namespace theoretica {

	/// @namespace theoretica::algebra
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

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.iat(i, j) = (i == j) ? nan() : 0;

			return m;
		}


		/// Overwrite the given vector with the error
		/// vector with NaN values. This function is used
		/// to signal an error.
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_error(Vector& v) {

			for (unsigned int i = 0; i < v.size(); ++i)
				v.iat(i) = nan();

			return v;
		}


		/// Overwrite a matrix with the identity matrix
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& mat_identity(Matrix& m) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.iat(i, j) = (i == j) ? 1 : 0;

			return m;
		}


		/// Overwrite a matrix with all zeroes
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& mat_zeroes(Matrix& m) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.iat(i, j) = 0;

			return m;
		}


		/// Overwrite a vector with all zeroes
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_zeroes(Vector& v) {

			for (unsigned int i = 0; i < v.size(); ++i)
				v.iat(i) = 0;

			return v;
		}


		/// Copy a matrix by overwriting another
		/// @param dest The matrix to overwrite
		/// @param src The matrix to copy
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& mat_copy(Matrix1& dest, Matrix2&& src) {

			if(src.rows() != dest.rows()) {
				TH_MATH_ERROR("algebra::mat_copy", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.cols()) {
				TH_MATH_ERROR("algebra::mat_copy", src.cols(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest.iat(i, j) = src.iget(i, j);

			return dest;
		}


		/// Copy a matrix by overwriting another.
		/// Equivalent to the operation dest = src
		/// @param dest The matrix to overwrite
		/// @param src The matrix to copy
		/// @return A reference to the overwritten matrix
		template<typename Vector1, typename Vector2>
		inline Vector1& vec_copy(Vector1& dest, Vector2&& src) {

			if(src.size() != dest.size()) {
				TH_MATH_ERROR("algebra::vec_copy", src.size(), INVALID_ARGUMENT);
				vec_error(dest);
				return dest;
			}

			for (unsigned int i = 0; i < src.size(); ++i)
				dest.iat(i) = src.iget(i);

			return dest;
		}


		/// Transpose the given matrix
		/// @param m The matrix to transpose
		/// @return A reference to the transposed matrix
		template<typename Matrix>
		inline Matrix& transpose(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::transpose", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			for (unsigned int i = 0; i < m.rows(); ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					
					auto buff = m.iat(i, j);
					m.iat(i, j) = m.iat(j, i);
					m.iat(j, i) = buff;
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
		inline Matrix1& transpose(Matrix1& dest, Matrix2&& src) {

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
					dest.iat(j, i) = src.iget(i, j);

			return dest;
		}


		/// Compute the hermitian of a given matrix
		/// and overwrite it.
		/// Equivalent to the operation m = m^H
		/// @param m The matrix to compute the hermitian of
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& hermitian(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::hermitian", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			for (unsigned int i = 0; i < m.rows(); ++i) {
				for (unsigned int j = 0; j < i; ++j) {
					
					auto buff = m.iat(i, j);
					m.iat(i, j) = conjugate(m.iat(j, i));
					m.iat(j, i) = conjugate(buff);
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
		inline Matrix1& hermitian(Matrix1& dest, Matrix2&& src) {

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
					dest.iat(j, i) = conjugate(src.iget(i, j));

			return dest;
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
			Matrix A, B;
			A.resize(m.rows(), m.cols());
			B.resize(m.rows(), m.cols());
			mat_copy(A, m);
			mat_identity(B);

			// Iterate on all columns
			for (unsigned int i = 0; i < m.rows(); ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < m.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < m.rows(); ++k) {
								A.iat(i, k) += A.iat(j, k);
								B.iat(i, k) += B.iat(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						TH_MATH_ERROR("algebra::invert", flag, IMPOSSIBLE_OPERATION);
						mat_error(m);
						return m;
					}
				}

				real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < m.rows(); ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.iat(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < m.rows(); ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
						B.iat(j, k) -= coeff * B.iat(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < m.rows(); ++j) {
					A.iat(i, j) *= inv_pivot;
					B.iat(i, j) *= inv_pivot;
				}
				
			}

			// Modify the matrix only when the inversion
			// has succeeded
			mat_copy(m, B);
			return m;
		}


		/// Invert the given matrix and overwrite it.
		/// Equivalent to the operation dest = src^-1
		/// @param dest The matrix to overwrite
		/// @param src The matrix to invert
		/// @return A reference to the inverted matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& inverse(Matrix1& dest, Matrix2&& src) {

			if(src.rows() != src.cols()) {
				TH_MATH_ERROR("algebra::inverse", src.rows(), INVALID_ARGUMENT);
				mat_error(dest);
				return dest;
			}

			// Prepare extended matrix (A|B)
			Matrix1 A;
			A.resize(src.rows(), src.cols());
			// dest.resize(src.rows(), src.cols());
			mat_copy(A, src);
			mat_identity(dest);

			// Iterate on all columns
			for (unsigned int i = 0; i < src.rows(); ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < src.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < src.rows(); ++k) {
								A.iat(i, k) += A.iat(j, k);
								dest.iat(i, k) += dest.iat(j, k);
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

				real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < src.rows(); ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.iat(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < src.rows(); ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
						dest.iat(j, k) -= coeff * dest.iat(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < src.rows(); ++j) {
					A.iat(i, j) *= inv_pivot;
					dest.iat(i, j) *= inv_pivot;
				}
				
			}

			return dest;
		}


		/// Compute the trace of the given matrix
		/// @param m A matrix of any type
		/// @return The trace of the matrix
		template<typename Matrix>
		inline auto trace(Matrix&& m) {

			auto sum = m.iget(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				sum += m.iget(i, i);

			return sum;
		}


		/// Compute the product of the elements of
		/// the main diagonal of a generic matrix
		/// @param m The input matrix
		/// @return The product of all the elements of the
		/// main diagonal of the input matrix
		template<typename Matrix>
		inline auto diagonal_product(Matrix&& m) {

			auto mul = m.iget(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				mul *= m.iget(i, i);

			return mul;
		}


		/// Compute the determinant of a matrix.
		/// The matrix must be square.
		/// Gauss Jordan elimination is used to reduce the
		/// matrix to a triangular matrix.
		/// @param m The matrix to compute the determinant of
		/// @return The determinant of the matrix
		template<typename Matrix>
		inline auto det(Matrix m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::det", m.rows(), INVALID_ARGUMENT);
				return nan();
			}

			typename std::remove_reference<Matrix>::type A;
			A.resize(m.rows(), m.cols());
			mat_copy(A, m);

			// Iterate on all columns
			for (unsigned int i = 0; i < m.rows(); ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.iat(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < m.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero.
						// The determinant does not change
						// when adding a row to another one
						if(A.iat(j, i) != 0) {

							for (unsigned int k = 0; k < m.rows(); ++k) {
								A.iat(i, k) += A.iat(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						return (decltype(m.iget(0, 0))) 0;
					}
				}

				const real inv_pivot = 1.0 / A.iat(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = i + 1; j < m.rows(); ++j) {

					// Multiplication coefficient for
					// the elision of Ajk
					const real coeff = A.iat(j, i) * inv_pivot;

					// The coefficient does not change
					// when adding a linear combination
					// of a row to another
					for (unsigned int k = 0; k < m.rows(); ++k) {
						A.iat(j, k) -= coeff * A.iat(i, k);
					}
				}
			}

			// The determinant of a (lower) triangular matrix
			// is the product of the elements on its diagonal
			return (decltype(m.iget(0, 0))) 1;
		}


		/// Return the determinant of a 2x2 matrix.
		/// @note No error checking is performed on the matrix size
		template<typename Matrix>
		inline real det_2x2(Matrix&& m) {
			return m.iget(0, 0) * m.iget(1, 1) - m.iget(1, 0) * m.iget(0, 1);
		}


		/// Return the determinant if the matrix is 3x3.
		/// @note No error checking is performed on the matrix size
		template<typename Matrix>
		inline real det_3x3(Matrix&& m) {
			return	m.iget(0, 0) * (m.iget(1, 1) * m.iget(2, 2) - m.iget(2, 1) * m.iget(1, 2)) -
					m.iget(0, 1) * (m.iget(1, 0) * m.iget(2, 2) - m.iget(2, 0) * m.iget(1, 2)) +
					m.iget(0, 2) * (m.iget(1, 0) * m.iget(2, 1) - m.iget(2, 0) * m.iget(1, 1));
		}


		/// Multiply a matrix by a scalar of any compatible type
		/// @param a A scalar value
		/// @param m The matrix to multiply
		/// @return A reference to the multiplied matrix
		template<typename Field, typename Matrix>
		inline Matrix& mat_scalmul(Field a, Matrix& m) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.iat(i, j) *= a;

			return m;
		}


		/// Multiply a matrix by a scalar of any compatible type
		/// which can be cast to the type of element of the output matrix.
		/// @param dest The matrix to overwrite with the result
		/// @param a A scalar value
		/// @param src The matrix to multiply
		/// @return A reference to the resulting matrix
		template<typename Field, typename Matrix1, typename Matrix2>
		inline Matrix1& mat_scalmul(Matrix1& dest, Field a, Matrix2&& src) {

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
					dest.iat(i, j) = a * src.iget(i, j);

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
		inline Vector& transform(Matrix&& A, Vector& v) {

			if(v.size() != A.cols()) {
				TH_MATH_ERROR("algebra::transform", v.size(), INVALID_ARGUMENT);
				vec_error(v);
				return v;
			}

			Vector res;
			res.resize(v.size());
			vec_zeroes(res);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res.iat(i) += A.iget(i, j) * v.iget(j);

			vec_copy(res, v);
			return v;
		}


		/// Apply a matrix transformation to a vector
		/// and store the result in the vector.
		/// Equivalent to the operation v = A * v
		/// @param res The matrix to overwrite with the result
		/// @param A The matrix transformation
		/// @param v The vector to transform
		/// @return A reference to the overwritten vector
		template<typename Matrix, typename Vector1, typename Vector2>
		inline Vector1& transform(Vector1& res, Matrix&& A, Vector2&& v) {

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
					res.iat(i) += A.iget(i, j) * v.iget(j);

			return res;
		}


		/// Construct a matrix with the given vector as diagonal
		/// and zeroes everywhere else
		/// @param The vector of diagonal elements
		/// @param res The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Vector, typename Matrix>
		inline Matrix& diagonal(Matrix& res, Vector&& v) {

			if(v.size() != res.cols()) {
				TH_MATH_ERROR("algebra::mat_diagonal", v.size(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			if(v.size() != res.rows()) {
				TH_MATH_ERROR("algebra::mat_diagonal", v.size(), INVALID_ARGUMENT);
				mat_error(res);
				return res;
			}

			for (unsigned int i = 0; i < res.rows(); ++i)
				for (unsigned int j = 0; j < res.cols(); ++j)
					res.iat(i, j) = (i == j) ? v.iget(i) : 0;

			return res;
		}


		// Operations involving multiple matrices


		/// Sum two matrices and store the result in the first matrix.
		/// Equivalent to the operation A = A + B
		/// @param A The first matrix to add and store the result
		/// @param B The second matrix to add
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix2& mat_sum(Matrix1& A, Matrix2&& B) {

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
					A.iat(i, j) = A.iget(i, j) + B.iget(i, j);

			return A;
		}


		/// Sum two matrices and store the result in another matrix
		/// Equivalent to the operation res = A + B
		/// @param A The first matrix to add
		/// @param B The second matrix to add
		/// @return A reference to the overwritten matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline Matrix1& mat_sum(Matrix1& res, Matrix2&& A, Matrix3&& B) {

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
					res.iat(i, j) = A.iget(i, j) + B.iget(i, j);

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
			Field1 alpha, Matrix1& A, Field2 beta, Matrix2&& B) {

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
					A.iat(i, j) = A.iget(i, j) * alpha + B.iget(i, j) * beta;

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
			Matrix1& res, Field1 alpha, Matrix2&& A, Field2 beta, Matrix3&& B) {

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
					res.iat(i, j) = A.iget(i, j) * alpha + B.iget(i, j) * beta;

			return A;
		}


		/// Multiply two matrices and store the result in the first matrix.
		/// Equivalent to the operation A = A * B
		/// @param A The first matrix to multiply and store the result
		/// @param B The second matrix to multiply
		/// @return A reference to the resulting matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& mat_mul(Matrix1& A, Matrix2&& B) {

			Matrix1 res;
			res.resize(A.rows(), B.cols());
			mat_zeroes(res);
			
			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < B.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						res.iat(i, j) += A.iget(i, k) * B.iget(k, j);

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
		inline Matrix1& mat_mul(Matrix1& res, Matrix2&& A, Matrix3&& B) {
			
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
						res.iat(i, j) += A.iget(i, k) * B.iget(k, j);

			return res;
		}


		// Linear transformations

		// LU decomposition

		// QR decomposition

		// Eigenvalues

		// Eigenvectors


	}

}

#endif
