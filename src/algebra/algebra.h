
///
/// @file algebra.h Linear algebra routines.
/// This file implements all linear algebra routines of the library,
/// working on templated data structures. The Matrix template must 
/// be a class with these methods:
/// - iget(i, j)  Get the element in the i-th row and j-th column
/// - iat(i, j)  Get a reference to the element in the i-th row and j-th column
/// - size()  Get the number of elements of matrix (N. rows * N. columns)
/// - rows()  Get the number of rows of the matrix
/// - cols()  Get the number of columns of the matrix
/// - resize() Change or set the size of the matrix
///

#ifndef THEORETICA_ALGEBRA_H
#define THEORETICA_ALGEBRA_H


namespace theoretica {


	/// Overwrite a matrix with the identity matrix
	template<typename Matrix>
	Matrix& mat_identity(Matrix& m) {

		for (unsigned int i = 0; i < m.rows(); ++i)
			for (unsigned int j = 0; j < m.cols(); ++j)
				m.iat(i, j) = (i == j) ? 1 : 0;

		return m;
	}


	/// Copy a matrix by overwriting another
	template<typename Matrix1, typename Matrix2>
	Matrix2& mat_copy(Matrix1&& src, Matrix2& dest) {

		if(src.rows() != dest.rows()) {
			TH_MATH_ERROR("th::mat_identity", src.rows(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		if(src.cols() != dest.cols()) {
			TH_MATH_ERROR("th::mat_identity", src.cols(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		for (unsigned int i = 0; i < src.rows(); ++i)
			for (unsigned int j = 0; j < src.cols(); ++j)
				dest.iat(i, j) = src.iat(i, j);

		return dest;
	}


	/// Overwrite the given matrix with the error
	/// matrix with NaN values on the diagonal and
	/// zeroes everywhere else.
	/// @param m The matrix to overwrite
	/// @return A reference to the overwritten matrix
	template<typename Matrix>
	Matrix& mat_error(Matrix& m) {

		for (unsigned int i = 0; i < m.rows(); ++i)
			for (unsigned int j = 0; j < m.cols(); ++j)
				m.iat(i, j) = (i == j) ? nan() : 0;

		return m;
	}


	/// Transpose the given matrix
	/// @param m The matrix to transpose
	/// @return A reference to the transposed matrix
	template<typename Matrix>
	Matrix& transpose(Matrix& m) {

		if(m.rows() != m.cols()) {
			TH_MATH_ERROR("th::transpose", m.rows(), INVALID_ARGUMENT);
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
	/// @param src The matrix to transpose
	/// @param dest The matrix to overwrite
	/// @return A reference to the overwritten matrix
	template<typename Matrix1, typename Matrix2>
	Matrix2& transposed(Matrix1&& src, Matrix2& dest) {

		// Check that the two matrices have the correct
		// number of rows and columns
		if(src.rows() != dest.cols()) {
			TH_MATH_ERROR("th::transposed", src.rows(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		if(src.cols() != dest.rows()) {
			TH_MATH_ERROR("th::transposed", src.cols(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		// Overwrite dest with the transpose of src
		for (unsigned int i = 0; i < src.rows(); ++i)
			for (unsigned int j = 0; j < src.cols(); ++j)
				dest.iat(j, i) = src.iat(i, j);

		return dest;
	}


	/// Invert the given matrix and overwrite it
	/// @param m The matrix to invert
	/// @return A reference to the inverted matrix
	template<typename Matrix>
	Matrix& invert(Matrix& m) {

		if(m.rows() != m.cols()) {
			TH_MATH_ERROR("th::invert", m.rows(), INVALID_ARGUMENT);
			mat_error(m);
			return m;
		}

		// Prepare extended matrix (A|B)
		Matrix A, B;
		A.resize(m.rows(), m.cols());
		B.resize(m.rows(), m.cols());
		mat_copy(m, A);
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
					TH_MATH_ERROR("th::invert", flag, IMPOSSIBLE_OPERATION);
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
		mat_copy(B, m);
		return m;
	}


	/// Invert the given matrix and overwrite it
	/// @param m The matrix to invert
	/// @return A reference to the inverted matrix
	template<typename Matrix1, typename Matrix2>
	Matrix2& inverse(Matrix1&& src, Matrix2& dest) {

		if(src.rows() != src.cols()) {
			TH_MATH_ERROR("th::invert", src.rows(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		// Prepare extended matrix (A|B)
		Matrix2 A;
		A.resize(src.rows(), src.cols());
		// dest.resize(src.rows(), src.cols());
		mat_copy(src, A);
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
					TH_MATH_ERROR("invert", flag, IMPOSSIBLE_OPERATION);
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


	/// Hermitian (conjugate transpose) of a matrix.
	/// The base type of the matrix needs to have
	/// a compatible conjugate() function.
	/// @param src The matrix to compute the hermitian of
	/// @param dest The matrix to overwrite
	/// @return A reference to the overwritten matrix
	template<typename Matrix1, typename Matrix2>
	Matrix2& hermitian(Matrix1&& src, Matrix2& dest) {

		// Check that the two matrices have the correct
		// number of rows and columns
		if(src.rows() != dest.cols()) {
			TH_MATH_ERROR("th::hermitian", src.rows(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		if(src.cols() != dest.rows()) {
			TH_MATH_ERROR("th::hermitian", src.cols(), INVALID_ARGUMENT);
			mat_error(dest);
			return dest;
		}

		// Overwrite dest with the transpose of src
		for (unsigned int i = 0; i < src.rows(); ++i)
			for (unsigned int j = 0; j < src.cols(); ++j)
				dest.iat(j, i) = conjugate(src.iat(i, j));

		return dest;
	}


	/// Compute the hermitian of a given matrix
	/// and overwrite it
	/// @param m The matrix to compute the hermitian of
	/// @return A reference to the overwritten matrix
	template<typename Matrix>
	Matrix& hermitian(Matrix& m) {

		if(m.rows() != m.cols()) {
			TH_MATH_ERROR("th::transpose", m.rows(), INVALID_ARGUMENT);
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


	/// Compute the trace of the given matrix
	/// @param m A matrix of any type
	/// @return The trace of the matrix
	template<typename Matrix>
	auto trace(Matrix&& m) {

		auto sum = m.iget(0, 0);
		const size_t n = min(m.rows(), m.cols());

		for (unsigned int i = 1; i < n; ++i)
			sum += m.iget(i, i);

		return sum;
	}


	template<typename Matrix>
	auto diagonal_product(Matrix&& m) {

		auto mul = m.iget(0, 0);
		const size_t n = min(m.rows(), m.cols());

		for (unsigned int i = 1; i < n; ++i)
			mul *= m.iget(i, i);

		return mul;
	}


	template<typename Matrix>
	auto det(Matrix&& m) {

		if(m.rows() != m.cols()) {
			TH_MATH_ERROR("th::det", m.rows(), INVALID_ARGUMENT);
			return nan();
		}

		Matrix A;
		A.resize(m.rows(), m.cols());
		mat_copy(m, A);

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

			real inv_pivot = 1.0 / A.iat(i, i);

			// Use the current row to make all other
			// elements of the column equal to zero
			for (unsigned int j = i + 1; j < m.rows(); ++j) {

				// Multiplication coefficient for
				// the elision of Ajk
				real coeff = A.iat(j, i) * inv_pivot;

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
		return (decltype(m.iget(0, 0))) diagonal_product(A);
	}


}

#endif
