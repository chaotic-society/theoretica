
///
/// @file algebra.h Linear algebra routines
///

#ifndef THEORETICA_ALGEBRA_H
#define THEORETICA_ALGEBRA_H


/// This file implements all linear algebra routines of the library,
/// working on templated data structures. The Matrix template must 
/// be a class with these methods:
/// - iget(i, j)  Get the element in the i-th row and j-th column
/// - iat(i, j)  Get a reference to the element in the i-th row and j-th column
/// - size()  Get the number of elements of matrix (N. rows * N. columns)
/// - rows()  Get the number of rows of the matrix
/// - cols()  Get the number of columns of the matrix


namespace theoretica {


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
	Matrix2& transposed(Matrix1& src, Matrix2& dest) {

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


	/// Hermitian (conjugate transpose) of a matrix.
	/// The base type of the matrix needs to have
	/// a compatible conjugate() function.
	/// @param src The matrix to compute the hermitian of
	/// @param dest The matrix to overwrite
	/// @return A reference to the overwritten matrix
	template<typename Matrix1, typename Matrix2>
	Matrix2& hermitian(Matrix1& src, Matrix2& dest) {

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


}

#endif
