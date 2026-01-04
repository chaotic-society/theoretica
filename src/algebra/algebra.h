
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

			using Type = vector_element_t<Vector>;

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

			using Type = vector_element_t<Vector>;

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
				TH_MATH_ERROR("algebra::mat_swap_rows", row1, MathError::InvalidArgument);
				return mat_error(A);
			}

			if (row2 >= A.rows()) {
				TH_MATH_ERROR("algebra::mat_swap_rows", row2, MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_swap_cols", col1, MathError::InvalidArgument);
				return mat_error(A);
			}

			if (col2 >= A.cols()) {
				TH_MATH_ERROR("algebra::mat_swap_cols", col2, MathError::InvalidArgument);
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


		/// Shift the diagonal of a matrix by the given amount, overwriting
		/// the matrix itself, as \f$(A + \sigma I)\f$.
		///
		/// @param A The matrix to shift the diagonal of
		/// @param sigma The amount to shift
		/// @return A reference to the modified matrix
		template<typename Matrix, typename Type = matrix_element_t<Matrix>>
		inline Matrix& mat_shift_diagonal(Matrix& A, const Type& sigma) {

			const unsigned int count = min(A.rows(), A.cols());

			for (unsigned int i = 0; i < count; ++i)
				A(i, i) += sigma;

			return A;
		}


		/// Compute the contribution of the inner product between
		/// a pair of elements of two vectors, automatically
		/// selecting whether to compute the conjugate or not.
		///
		/// @param v_i The first element of the pair
		/// @param w_i The second element of the pair, which will be conjugated if needed
		/// @return The contribution of the element pair
		template<typename Type>
		inline auto pair_inner_product(const Type& v_i, const Type& w_i) {
			return v_i * w_i;
		}

		/// Compute the contribution of the inner product between
		/// a pair of elements of two vectors, automatically
		/// selecting whether to compute the conjugate or not.
		///
		/// @param v_i The first element of the pair
		/// @param w_i The second element of the pair, which will be conjugated if needed
		/// @return The contribution of the element pair
		template<typename Type>
		inline auto pair_inner_product(const complex<Type>& v_i, const complex<Type>& w_i) {
			return v_i * conjugate(w_i);
		}


		/// Returns the square of the Euclidean/Hermitian norm
		/// of the given vector
		///
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto sqr_norm(const Vector& v) {

			auto sum = pair_inner_product(v[0], v[0]);

			// Use conjugation for complex numbers
			for (unsigned int i = 1; i < v.size(); ++i)
				sum += pair_inner_product(v[i], v[i]);

			return sum;
		}


		/// Returns the Euclidean/Hermitian norm of the given vector
		///
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto norm(const Vector& v) {

			return vector_element_t<Vector>(sqrt(sqr_norm(v)));
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
				TH_MATH_ERROR("algebra::normalize", m, MathError::DivByZero);
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
				TH_MATH_ERROR("algebra::make_normalized", m, MathError::DivByZero);
				vec_error(v);
				return v;
			}

			for (unsigned int i = 0; i < v.size(); ++i)
				v[i] /= m;

			return v;
		}


		/// Computes the dot product between two vectors,
		/// using the Hermitian form if needed.
		///
		/// @param v The first vector
		/// @param w The second vector
		/// @return The dot product of the two vectors
		template<typename Vector1, typename Vector2>
		inline auto dot(const Vector1& v, const Vector2& w) {

			if(v.size() != w.size()) {
				TH_MATH_ERROR("algebra::dot", v.size(), MathError::InvalidArgument);
				return vector_element_t<Vector1>(nan());
			}

			auto sum = pair_inner_product(v[0], w[0]);

			// Use conjugation for complex numbers
			for (unsigned int i = 1; i < v.size(); ++i)
				sum += pair_inner_product(v[i], w[i]);

			return sum;
		}


		/// Compute the cross product between two tridimensional vectors
		/// @param v1 The first tridimensional vector
		/// @param w The second tridimensional vector
		/// @return The cross product of the two vectors
		template<typename Vector1, typename Vector2>
		inline Vector1 cross(const Vector1& v1, const Vector2& v2) {

			Vector1 v3;
			v3.resize(3);

			if(v1.size() != 3) {
				TH_MATH_ERROR("algebra::cross", v1.size(), MathError::InvalidArgument);
				vec_error(v3);
				return v3;
			}

			if(v2.size() != 3) {
				TH_MATH_ERROR("algebra::cross", v2.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::make_transposed", m.rows(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::transpose", src.rows(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.rows()) {
				TH_MATH_ERROR("algebra::transpose", src.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::hermitian", m.rows(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::hermitian", src.rows(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.rows()) {
				TH_MATH_ERROR("algebra::hermitian", src.cols(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			// Overwrite dest with the transpose of src
			for (unsigned int i = 0; i < src.rows(); ++i)
				for (unsigned int j = 0; j < src.cols(); ++j)
					dest(j, i) = conjugate(src(i, j));

			return dest;
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
				TH_MATH_ERROR("algebra::mat_scalmul", src.rows(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			if(src.cols() != dest.cols()) {
				TH_MATH_ERROR("algebra::mat_scalmul", src.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::vec_scalmul", src.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::apply_transform", v.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::transform", v.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::transform", v.size(), MathError::InvalidArgument);
				vec_error(res);
				return res;
			}

			if(res.size() != v.size()) {
				TH_MATH_ERROR("algebra::transform", res.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_sum", A.rows(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", A.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_sum", A.rows(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", A.cols(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_sum", res.rows(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_sum", res.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_diff", A.rows(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", A.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_diff", A.rows(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", A.cols(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_diff", res.rows(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_diff", res.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_lincomb", A.rows(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.cols(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::mat_lincomb", A.rows(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", A.cols(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			if(res.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_lincomb", res.rows(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			if(res.cols() != A.cols()) {
				TH_MATH_ERROR("algebra::mat_lincomb", res.cols(), MathError::InvalidArgument);
				mat_error(res);
				return res;
			}

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					res(i, j) = A(i, j) * alpha + B(i, j) * beta;

			return A;
		}


		/// Multiply two matrices and store the result in the first matrix,
		/// equivalent to the operation \f$R = A B\f$.
		///
		/// @param A The first matrix to multiply and store the result
		/// @param B The second matrix to multiply
		/// @return The resulting matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3 = Matrix1>
		inline Matrix3 mat_mul(const Matrix1& A, const Matrix2& B) {

			Matrix3 R;
			R.resize(A.rows(), B.cols());

			if(A.cols() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", A.cols(), MathError::InvalidArgument);
				mat_error(R);
				return R;
			}

			mat_zeroes(R);
			
			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < B.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						R(i, j) += A(i, k) * B(k, j);

			return R;
		}


		/// Multiply two matrices and store the result in another matrix,
		/// equivalent to the operation \f$R = A B\f$.
		///
		/// @param R The matrix to overwrite with the result
		/// @param A The first matrix to multiply
		/// @param B The second matrix to multiply
		/// @return A reference to the modified matrix
		template<typename Matrix1, typename Matrix2, typename Matrix3>
		inline Matrix1& mat_mul(Matrix1& R, const Matrix2& A, const Matrix3& B) {
			
			if(R.rows() != A.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", R.rows(), MathError::InvalidArgument);
				mat_error(R);
				return R;
			}

			if(R.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_mul", R.cols(), MathError::InvalidArgument);
				mat_error(R);
				return R;
			}

			if(A.cols() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_mul", A.cols(), MathError::InvalidArgument);
				mat_error(R);
				return R;
			}

			mat_zeroes(R);

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < B.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						R(i, j) += A(i, k) * B(k, j);

			return R;
		}


		/// Multiply the transpose of a matrix by another matrix,
		/// equivalent to the operation \f$R = A^T B\f$.
		///
		/// @note This function is faster then writing algebra::transpose(A) * B
		/// and should be preferred.
		///
		/// @param A The matrix to transpose and then multiply
		/// @param B The second matrix to multiply by
		/// @return The result of the multiplication by the transpose of the first
		/// matrix and the second matrix.
		template<typename Matrix1, typename Matrix2, typename Matrix3 = Matrix1>
		inline Matrix3 mat_transpose_mul(const Matrix1& A, const Matrix2& B) {

			Matrix3 R;
			R.resize(A.cols(), B.cols());

			if(A.rows() != B.rows()) {
				TH_MATH_ERROR("algebra::mat_transpose_mul", A.rows(), MathError::InvalidArgument);
				return mat_error(R);
			}

			mat_zeroes(R);

			for (unsigned int i = 0; i < R.rows(); ++i)
				for (unsigned int j = 0; j < R.cols(); ++j)
					for (unsigned int k = 0; k < A.rows(); ++k)
						R(i, j) += A(k, i) * B(k, j);

			return R;
		}


		/// Multiply a matrix by the transpose of another matrix,
		/// equivalent to the operation \f$R = A B^T\f$.
		///
		/// @note This function is faster then writing A * algebra::transpose(B)
		/// and should be preferred.
		///
		/// @param A The first matrix to multiply
		/// @param B The second matrix to transpose and multiply by
		/// @return The result of the multiplication by the first
		/// matrix and the transpose of second matrix.
		template<typename Matrix1, typename Matrix2, typename Matrix3 = Matrix1>
		inline Matrix3 mat_mul_transpose(const Matrix1& A, const Matrix2& B) {

			Matrix3 R;
			R.resize(A.rows(), B.rows());

			if(A.cols() != B.cols()) {
				TH_MATH_ERROR("algebra::mat_mul_transpose", A.rows(), MathError::InvalidArgument);
				return mat_error(R);
			}

			mat_zeroes(R);

			for (unsigned int i = 0; i < R.rows(); ++i)
				for (unsigned int j = 0; j < R.cols(); ++j)
					for (unsigned int k = 0; k < A.cols(); ++k)
						R(i, j) += A(i, k) * B(j, k);

			return R;
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

			for (unsigned int i = 0; i < A.rows(); ++i) {
				for (unsigned int j = 0; j < A.cols(); ++j) {
			
					const auto diff = abs(A(i, j) - B(i, j));

					if(diff > tolerance || is_nan(diff))
						return false;
				}
			}

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
				TH_MATH_ERROR("algebra::vec_sum", v1.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::vec_sum", v1.size(), MathError::InvalidArgument);
				vec_error(res);
				return res;
			}

			if(res.size() != v1.size()) {
				TH_MATH_ERROR("algebra::vec_sum", res.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::vec_diff", v1.size(), MathError::InvalidArgument);
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
				TH_MATH_ERROR("algebra::vec_diff", v1.size(), MathError::InvalidArgument);
				vec_error(res);
				return res;
			}

			if(res.size() != v1.size()) {
				TH_MATH_ERROR("algebra::vec_diff", res.size(), MathError::InvalidArgument);
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
		/// in the comparison, defaults to ALGEBRA_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_symmetric(const Matrix& m, real tolerance = ALGEBRA_ELEMENT_TOL) {

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
		/// in the comparison, defaults to ALGEBRA_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_lower_triangular(const Matrix& m, real tolerance = ALGEBRA_ELEMENT_TOL) {

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
		/// in the comparison, defaults to ALGEBRA_ELEMENT_TOL
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_upper_triangular(const Matrix& m, real tolerance = ALGEBRA_ELEMENT_TOL) {

			if(!is_square(m))
				return false;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < i; ++j)
					if (abs(m(i, j)) > tolerance)
						return false;

			return true;
		}


		/// Compute the density of a matrix, counting the proportion
		/// of non-zero (bigger in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param A The matrix to compute the density of
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of non-zero elements of the matrix.
		template<typename Matrix, enable_matrix<Matrix> = true>
		inline real density(const Matrix& A, real tolerance = 1E-12) {

			long unsigned int n = 0;

			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = 0; j < A.cols(); ++j)
					if (abs(A(i, j)) > tolerance)
						n++;

			return (real(n) / A.rows()) / A.cols();
		}


		/// Compute the sparsity of a matrix, counting the proportion
		/// of zero (smaller in module than the given tolerance) elements
		/// with respect to the total number of elements.
		///
		/// @param A The matrix to compute the sparsity of
		/// @param tolerance The minimum tolerance in absolute value
		/// to consider an element non-zero.
		/// @return A real number between 0 and 1 representing the
		/// proportion of zero elements of the matrix.
		template<typename Matrix, enable_matrix<Matrix> = true>
		inline real sparsity(const Matrix& A, real tolerance = 1E-12) {

			return 1.0 - density(A, tolerance);
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

			// Check the shapes of A, L and U
			unsigned int err = 0;

			// Make sure to allocate L and U if they are empty
			L.resize(A.rows(), A.cols());
			U.resize(A.rows(), A.cols());

			if (!is_square(A))
				err = A.rows();

			if (A.rows() != L.rows())
				err = L.rows();

			if (A.cols() != L.cols())
				err = L.cols();

			if (A.rows() != U.rows())
				err = U.rows();

			if (A.cols() != U.cols())
				err = U.cols();

			if (err) {
				TH_MATH_ERROR("algebra::decompose_lu", err, MathError::InvalidArgument);
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
						U(i, j) -= pair_inner_product(L(i, k), U(k, j));
				}

				for(unsigned int j = i + 1; j < A.rows(); ++j) {

					L(j, i) = A(j, i);

					for(unsigned int k = 0; k < i; ++k)
						L(j, i) -= pair_inner_product(L(j, k), U(k, i));

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
				TH_MATH_ERROR("algebra::decompose_lu_inplace", A.rows(), MathError::InvalidArgument);
				mat_error(A);
				return A;
			}

			for(unsigned int j = 0; j < A.cols(); ++j) {

				for(unsigned int i = j + 1; i < A.rows(); ++i) {

					A(i, j) /= A(j, j);

					for(unsigned int k = j + 1; k < A.rows(); ++k)
						A(i, k) -= pair_inner_product(A(i, j), A(j, k));
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
				TH_MATH_ERROR("algebra::decompose_cholesky", A.rows(), MathError::InvalidArgument);
				return mat_error(L);
			}

			if (!is_symmetric(A)) {
				TH_MATH_ERROR("algebra::decompose_cholesky", false, MathError::InvalidArgument);
				return mat_error(L);
			}

			mat_zeroes(L);

			for (unsigned int i = 0; i < A.cols(); ++i) {
				
				for (unsigned int j = 0; j <= i; ++j) {
					
					Type sum = pair_inner_product(L(i, 0), L(j, 0));

					for (unsigned int k = 1; k < j; ++k)
						sum += pair_inner_product(L(i, k), L(j, k));

					if (i == j) {

						const Type sqr_diag = A(j, j) - sum;

						// Additional check to ensure that the matrix is positive definite
						if (sqr_diag < MACH_EPSILON) {
							TH_MATH_ERROR("algebra::decompose_cholesky", sqr_diag, MathError::InvalidArgument);
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


		/// Decompose a symmetric positive definite matrix in-place,
		/// overwriting the starting matrix, without using additional space.
		///
		/// @param A The symmetric, positive definite matrix to decompose and overwrite
		/// with the lower triangular matrix.
		template<typename Matrix>
		inline Matrix& decompose_cholesky_inplace(Matrix& A) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("algebra::decompose_cholesky_inplace", A.rows(), MathError::InvalidArgument);
				return mat_error(A);
			}

			if (!is_symmetric(A)) {
				TH_MATH_ERROR("algebra::decompose_cholesky_inplace", false, MathError::InvalidArgument);
				return mat_error(A);
			}

			// Compute the Cholesky decomposition in-place
			for (unsigned int k = 0; k < A.rows(); ++k) {

				A(k, k) = sqrt(A(k, k));

				for (unsigned int i = k + 1; i < A.rows(); ++i)
					A(i, k) = A(i, k) / A(k, k);

				for (unsigned int j = k + 1; j < A.rows(); ++j)
					for (unsigned int i = j; i < A.cols(); ++i)
						A(i, j) -= pair_inner_product(A(i, k), A(j, k));
			}

			// Zero out elements over the diagonal
			for (unsigned int i = 0; i < A.rows(); ++i)
				for (unsigned int j = i + 1; j < A.rows(); ++j)
					A(i, j) = Type(0.0);
					

			return A;
		}


		// Linear system solvers


		/// Solve the linear system \f$L \vec x = b\f$ for lower triangular \f$L\f$.
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
				TH_MATH_ERROR("algebra::solve_triangular_lower", false, MathError::InvalidArgument);
				return vec_error(x);
			}

			if (b.size() != L.rows()) {
				TH_MATH_ERROR("algebra::solve_triangular_lower", b.size(), MathError::InvalidArgument);
				return vec_error(x);
			}

			// Solve using forward substitution
			for (unsigned int i = 0; i < L.cols(); ++i) {
				
				Type sum = L(i, 0) * x[0];

				for (unsigned int j = 1; j < i; ++j)
					sum += L(i, j) * x[j];

				if (abs(L(i, i)) < MACH_EPSILON) {
					TH_MATH_ERROR("algebra::solve_triangular_lower", L(i, i), MathError::DivByZero);
					return vec_error(x);
				}

				x[i] = (b[i] - sum) / L(i, i);
			}

			return x;
		}


		/// Solve the linear system \f$U \vec x = b\f$ for upper triangular \f$U\f$.
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
				TH_MATH_ERROR("solve_triangular_upper", false, MathError::InvalidArgument);
				return vec_error(x);
			}

			if (b.size() != U.rows()) {
				TH_MATH_ERROR("solve_triangular_upper", b.size(), MathError::InvalidArgument);
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
					TH_MATH_ERROR("solve_triangular_upper", U(i, i), MathError::DivByZero);
					return vec_error(x);
				}

				x[i] = (b[i] - sum) / U(i, i);
			}

			return x;
		}


		/// Solve the linear system \f$T \vec x = b\f$ for triangular \f$T\f$.
		/// The correct solver is selected depending on the elements of \f$T\f$,
		/// if the property of the matrix is known a priori, calling the
		/// specific function is more efficient.
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


		/// Solve the linear system \f$A \vec x = \vec b\f$, finding \f$\vec x\f$,
		/// where the matrix A has **already undergone in-place LU decomposition**.
		/// Forward and backward elimination is used to solve the system in place.
		/// This routine is particularly efficient for solving linear systems
		/// multiple times with the same matrix but different vectors.
		/// The input vector is overwritten, as to not use any additional memory.
		/// 
		/// @param A The matrix of the linear system, after in-place LU decomposition
		/// @param b The known vector, to be overwritten with the solution
		/// @return A reference to the overwritten vector solution
		template<typename Matrix, typename Vector>
		inline Vector& solve_lu_inplace(const Matrix& A, Vector& b) {

			if (!is_square(A)) {
				TH_MATH_ERROR("algebra::solve_lu_inplace", A.rows(), MathError::InvalidArgument);
				return vec_error(b);
			}

			if (A.rows() != b.size()) {
				TH_MATH_ERROR("algebra::solve_lu_inplace", A.rows(), MathError::InvalidArgument);
				return vec_error(b);
			}

			using Type = matrix_element_t<Matrix>;

			// Forward elimination: solves the lower system (L matrix).
			for (unsigned int i = 1; i < A.rows(); ++i) {
				
				Type sum = Type(0.0);

				for (unsigned int j = 0; j < i; ++j)
					sum += A(i, j) * b[j];

				b[i] -= sum;
			}

			// Backward elimination: solves the upper system (U matrix).
			for (int i = A.rows() - 1; i >= 0; --i) {

				Type sum = Type(0.0);
				
				for (unsigned int j = i + 1; j < A.rows(); ++j)
					sum += A(i, j) * b[j];

				b[i] = (b[i] - sum) / A(i, i);
			}
		
			return b;
		}


		/// Solve the linear system \f$A \vec x = \vec b\f$, finding \f$\vec x\f$.
		/// In-place LU decomposition is used on A, followed by forward and backward
		/// elimination.
		/// 
		/// @param A The matrix of the linear system
		/// @param b The known vector
		/// @return The unknown vector
		template<typename Matrix, typename Vector>
		inline Vector solve_lu(Matrix A, Vector b) {

			// Apply in-place LU decomposition
			decompose_lu_inplace(A);

			// Apply forward and backward substitution
			return solve_lu_inplace(A, b);
		}


		/// Use the LU decomposition of a matrix to solve its associated linear system,
		/// solving \f$A \vec x = \vec b\f$ for \f$\vec b\f$. When solving the same linear
		/// system over and over, it is advantageous to compute its LU decomposition
		/// using decompose_lu and then use the decomposition to solve the system for
		/// different known vectors, reducing the overall computational cost.
		///
		/// @param L The lower triangular matrix
		/// @param U The upper triangular matrix
		/// @param b The known vector
		/// @return The vector solution \f$\vec x\f$.
		template<typename Matrix1, typename Matrix2, typename Vector>
		inline Vector solve_lu(const Matrix1& L, const Matrix2& U, const Vector& b) {

			Vector x = b;

			if (!is_square(L)) {
				TH_MATH_ERROR("algebra::solve_lu", L.rows(), MathError::InvalidArgument);
				return vec_error(x);
			}

			if (!is_square(U)) {
				TH_MATH_ERROR("algebra::solve_lu", U.rows(), MathError::InvalidArgument);
				return vec_error(x);
			}

			if (L.rows() != U.rows()) {
				TH_MATH_ERROR("algebra::solve_lu", U.rows(), MathError::InvalidArgument);
				return vec_error(x);
			}

			if (b.size() != L.rows()) {
				TH_MATH_ERROR("algebra::solve_lu", b.size(), MathError::InvalidArgument);
				return vec_error(x);
			}

			using Type1 = matrix_element_t<Matrix1>;

			// Forward elimination for L
			for (unsigned int i = 1; i < L.rows(); ++i) {
				
				Type1 sum = Type1(0.0);

				for (unsigned int j = 0; j < i; ++j)
					sum += L(i, j) * x[j];

				x[i] -= sum;
			}

			using Type2 = matrix_element_t<Matrix2>;

			// Backward elimination for U
			for (int i = U.rows() - 1; i >= 0; --i) {

				Type2 sum = Type2(0.0);
				
				for (unsigned int j = i + 1; j < U.rows(); ++j)
					sum += U(i, j) * x[j];

				if (abs(U(i, i)) < MACH_EPSILON) {
					TH_MATH_ERROR("algebra::solve_lu", U(i, i), MathError::DivByZero);
					return vec_error(x);
				}

				x[i] = (x[i] - sum) / U(i, i);
			}

			return x;
		}


		/// Solve a linear system \f$A \vec x = \vec b\f$ defined by a symmetric
		/// positive definite matrix, using the Cholesky decomposition \f$L\f$ constructed
		/// so that \f$A = LL^T\f$.
		///
		/// @param L The already computed Cholesky decomposition of the symmetric
		/// positive definite matrix describing the system.
		/// @param b The known vector
		/// @return The unknown vector \f$\vec x\f$
		template<typename Matrix, typename Vector>
		inline Vector solve_cholesky(const Matrix& L, const Vector& b) {

			Vector x = b;

			if (!is_square(L)) {
				TH_MATH_ERROR("algebra::solve_cholesky", L.rows(), MathError::InvalidArgument);
				return vec_error(x);
			}

			if (L.rows() != b.size()) {
				TH_MATH_ERROR("algebra::solve_cholesky", b.size(), MathError::InvalidArgument);
				return vec_error(x);
			}

			using Type = matrix_element_t<Matrix>;

			// Forward elimination for L
			for (unsigned int i = 0; i < L.rows(); ++i) {
				
				Type sum = Type(0.0);

				for (unsigned int j = 0; j < i; ++j)
					sum += L(i, j) * x[j];

				x[i] = (x[i] - sum) / L(i, i);
			}

			// Backward elimination for L transpose
			for (int i = L.rows() - 1; i >= 0; --i) {

				Type sum = Type(0.0);
				
				for (unsigned int j = i + 1; j < L.rows(); ++j)
					sum += L(j, i) * x[j];

				x[i] = (x[i] - sum) / L(i, i);
			}

			return x;
		}


		/// Solve the linear system \f$A \vec x = \vec b\f$, finding \f$\vec x\f$
		/// using the best available algorithm.
		/// 
		/// @param A The matrix of the linear system
		/// @param b The known vector
		/// @return The unknown vector \f$\vec x\f$
		template<typename Matrix, typename Vector>
		inline Vector solve(const Matrix& A, const Vector& b) {

			// Use LU decomposition to solve the linear system
			return solve_lu(A, b);
		}


		// Other composite operations


		/// Invert the given matrix.
		/// Equivalent to the operation dest = src^-1
		/// @param dest The matrix to overwrite
		/// @param src The matrix to invert
		/// @return A reference to the inverted matrix
		template<typename Matrix1, typename Matrix2>
		inline Matrix1& inverse(Matrix1& dest, const Matrix2& src) {

			if(src.rows() != src.cols()) {
				TH_MATH_ERROR("algebra::inverse", src.rows(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			if(dest.rows() != src.rows()) {
				TH_MATH_ERROR("algebra::inverse", dest.rows(), MathError::InvalidArgument);
				mat_error(dest);
				return dest;
			}

			if(dest.cols() != src.cols()) {
				TH_MATH_ERROR("algebra::inverse", dest.cols(), MathError::InvalidArgument);
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
						TH_MATH_ERROR("algebra::inverse", flag, MathError::ImpossibleOperation);
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
		/// Equivalent to the operation \f$m^-1\f$
		///
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
		///
		/// @param m The matrix to invert
		/// @return A reference to the inverted matrix
		template<typename Matrix>
		inline Matrix& invert(Matrix& m) {

			if(m.rows() != m.cols()) {
				TH_MATH_ERROR("algebra::invert", m.rows(), MathError::InvalidArgument);
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


		/// Compute the determinant of a square matrix.
		/// In-place LU decomposition is used to reduce the
		/// matrix to triangular form.
		///
		/// @param A The matrix to compute the determinant of
		/// @return The determinant of the matrix
		template<typename Matrix>
		inline auto det(const Matrix& A) {

			Matrix LU = A;
			decompose_lu_inplace(LU);

			// The determinant of a triangular matrix
			// is the product of the elements on its diagonal
			return diagonal_product(LU);
		}


		// Eigensolvers


		/// Compute the Rayleigh quotient \f$\frac{x^T A x}{x^T x}\f$ of a vector
		/// with respect to a matrix. This value is particularly
		/// useful in the context of eigensolvers.
		///
		/// @param A The matrix
		/// @param x The vector
		/// @return The Rayleigh quotient of x with respect to A
		template<typename Matrix, typename Vector>
		inline auto rayleigh_quotient(const Matrix& A, const Vector& x) {

			const auto p = dot(x, x);

			// Check for division by zero
			if (abs(p) < MACH_EPSILON) {
				TH_MATH_ERROR("rayleigh_quotient", abs(p), MathError::DivByZero);
				return vector_element_t<Vector>(nan());
			}

			return dot(x, transform(A, x)) / p;
		}


		/// Find the biggest eigenvalue in module \f$|\lambda_i|\f$ of a
		/// square matrix using the power method (Von Mises iteration).
		///
		/// @param A The matrix to find the biggest eigenvalue of
		/// @param x The starting vector (a random vector is a good choice)
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return The biggest eigenvalue of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector>
		inline auto eigenvalue_power(
			const Matrix& A, const Vector& x,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenvalue_power", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (x.size() != A.rows()) {
				TH_MATH_ERROR("eigenvalue_power", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Apply a first iteration to initialize
			// the current and previous vectors
			Vector x_prev = x;
			Vector x_curr = normalize(A * x_prev);

			// Iteration counter
			unsigned int i;

			// Iteratively apply the matrix to the vector
			for (i = 1; i <= max_iter; ++i) {
				
				x_prev = x_curr;
				x_curr = normalize(transform(A, x_prev));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			// The algorithm did not converge
			if (i > max_iter) {
				TH_MATH_ERROR("eigenvalue_power", i, MathError::NoConvergence);
				return Type(nan());
			}

			return dot(x_curr, A * x_curr);
		}


		/// Find the biggest eigenvalue in module \f$|\lambda_i|\f$ of a
		/// square matrix and its corresponding eigenvector (eigenpair),
		/// using the power method (Von Mises iteration).
		///
		/// @param A The matrix to find the biggest eigenvalue of
		/// @param x The starting vector (a random vector is a good choice)
		/// @param v The vector to overwrite with the eigenvector
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return The biggest eigenvalue of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector1, typename Vector2 = Vector1>
		inline auto eigenpair_power(
			const Matrix& A, const Vector1& x, Vector2& v,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenpair_power", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (x.size() != A.rows()) {
				TH_MATH_ERROR("eigenpair_power", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			if (v.size() != x.size()) {
				TH_MATH_ERROR("eigenpair_power", v.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Apply a first iteration to initialize
			// the current and previous vectors
			Vector1 x_prev = x;
			Vector1 x_curr = normalize(A * x_prev);

			// Iteration counter
			unsigned int i;

			// Iteratively apply the matrix to the vector
			for (i = 1; i <= max_iter; ++i) {
				
				x_prev = x_curr;
				x_curr = normalize(transform(A, x_prev));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			// The algorithm did not converge
			if (i > max_iter) {
				TH_MATH_ERROR("eigenpair_power", i, MathError::NoConvergence);
				return Type(nan());
			}

			// Overwrite with the eigenvector
			vec_copy(v, x_curr);

			return dot(x_curr, A * x_curr);
		}


		/// Find the eigenvalue with the smallest inverse \f$\frac{1}{|\lambda_i|}\f$ of
		/// a square matrix, using the inverse power method with parameter equal to 0.
		///
		/// @param A The matrix to find the smallest eigenvalue of
		/// @param x The starting vector (a random vector is a good choice)
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return The eigenvalue with the smallest inverse of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector>
		inline auto eigenvalue_inverse(
			const Matrix& A, const Vector& x,
			real tolerance = ALGEBRA_EIGEN_TOL,
			unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenvalue_inverse", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (A.rows() != x.size()) {
				TH_MATH_ERROR("eigenvalue_inverse", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Compute the LU decomposition of A to speed up system solution
			Matrix LU = A;
			decompose_lu_inplace(LU);

			// Compute the first step to initialize the two vectors
			Vector x_prev = normalize(x);
			Vector x_curr = x_prev;
			solve_lu_inplace(LU, x_curr);

			// Iteration counter
			unsigned int i;

			for (i = 1; i <= max_iter; ++i) {
				
				// Solve the linear system using the LU decomposition
				x_prev = x_curr;
				make_normalized(solve_lu_inplace(LU, x_curr));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			// The algorithm did not converge
			if (i > max_iter) {
				TH_MATH_ERROR("eigenvalue_inverse", i, MathError::NoConvergence);
				return Type(nan());
			}

			// A and inverse(A) have the same eigenvectors
			return dot(x_curr, A * x_curr);
		}


		/// Find the eigenvalue with the smallest inverse \f$\frac{1}{|\lambda_i|}\f$
		/// of a square matrix and its corresponding eigenvector, using the
		/// inverse power method with parameter equal to 0.
		///
		/// @param A The matrix to find the smallest eigenvalue of
		/// @param x The starting vector (a random vector is a good choice)
		/// @param v The vector to overwrite with the eigenvector
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return The eigenvalue with the smallest inverse of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector1, typename Vector2>
		inline auto eigenpair_inverse(
			const Matrix& A, const Vector1& x, Vector2& v,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenpair_inverse", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (A.rows() != x.size()) {
				TH_MATH_ERROR("eigenpair_inverse", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Compute the LU decomposition of A to speed up system solution
			Matrix LU = A;
			decompose_lu_inplace(LU);

			// Compute the first step to initialize the two vectors
			Vector1 x_prev = normalize(x);
			Vector1 x_curr = x_prev;
			solve_lu_inplace(LU, x_curr);

			// Iteration counter
			unsigned int i;

			for (i = 1; i <= max_iter; ++i) {
				
				// Solve the linear system using the LU decomposition
				x_prev = x_curr;
				make_normalized(solve_lu_inplace(LU, x_curr));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			// The algorithm did not converge
			if (i > max_iter) {
				TH_MATH_ERROR("eigenpair_inverse", i, MathError::NoConvergence);
				return Type(nan());
			}

			// Overwrite with the eigenvector
			vec_copy(v, x_curr);

			// A and inverse(A) have the same eigenvectors
			return dot(x_curr, A * x_curr);
		}


		/// Find the eigenvector associated with a given
		/// approximated eigenvalue using the inverse power method.
		///
		/// @note The algorithm is unstable when the approximation
		/// of the eigenvalue is too close to the true value. A value
		/// not too close to the actual eigenvalue and far away from
		/// the other eigenvalues should be used.
		///
		/// @param A The matrix to find the eigenvector of
		/// @param lambda The eigenvalue
		/// @param x The starting vector (a random vector is a good choice)
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return The approximate eigenvector associated with the eigenvalue
		template<typename Matrix, typename Vector, typename T = matrix_element_t<Matrix>>
		inline Vector eigenvector_inverse(
			const Matrix& A, const T& lambda, const Vector& x,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenvector_inverse", is_square(A), MathError::InvalidArgument);
				Vector x;
				return vec_error(x);
			}

			if (A.rows() != x.size()) {
				TH_MATH_ERROR("eigenvector_inverse", x.size(), MathError::InvalidArgument);
				Vector x;
				return vec_error(x);
			}

			Matrix LU = A;
			mat_shift_diagonal(LU, -lambda);

			// Compute the LU decomposition of A
			// to speed up system solution
			decompose_lu_inplace(LU);

			// Compute the first step to initialize the two vectors
			Vector v_prev = normalize(x);
			Vector v_curr = v_prev;
			solve_lu_inplace(LU, v_curr);

			// Iteration counter
			unsigned int i;

			for (i = 1; i <= max_iter; ++i) {
				
				// Solve the linear system using the LU decomposition
				v_prev = v_curr;
				make_normalized(solve_lu_inplace(LU, v_curr));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(v_curr - v_prev) <= tolerance ||
					norm(v_curr + v_prev) <= tolerance)
					break;
			}

			// The algorithm did not converge
			if (i > max_iter) {
				TH_MATH_ERROR("eigenvector_inverse", i, MathError::NoConvergence);
				return vec_error(v_curr);
			}

			return v_curr;
		}


		/// Compute an eigenvalue of a square matrix using
		/// the Rayleigh quotient iteration method.
		///
		/// @param A The matrix to compute the eigenvalue of
		/// @param lambda The starting value for the Rayleigh quotient
		/// @param x The starting approximation for the eigenvector
		/// (a random vector is a good choice).
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return An approximate eigenvalue of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector, typename T = matrix_element_t<Matrix>>
		inline auto eigenvalue_rayleigh(
			const Matrix& A, const T& lambda, const Vector& x,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenvalue_rayleigh", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (A.rows() != x.size()) {
				TH_MATH_ERROR("eigenvalue_rayleigh", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Keep track of the shifted matrix
			Matrix A_shift = A;
			mat_shift_diagonal(A_shift, -lambda);

			Type lambda_prev = lambda;
			Type lambda_curr = lambda;

			Vector x_prev = normalize(x);
			Vector x_curr = normalize(solve(A_shift, x));

			unsigned int i;

			for (i = 1; i <= max_iter; ++i) {
					
				// Update the eigenvalue approximation
				// using Rayleigh quotient (avoiding normalization)
				lambda_prev = lambda_curr;
				lambda_curr = dot(x_curr, transform(A, x_curr));

				// Shift the diagonal by the difference between
				// subsequent eigenvalues steps, to avoid copying matrix A
				mat_shift_diagonal(A_shift, lambda_prev - lambda_curr);

				// Solve the linear system each time, as the shifted
				// matrix is continuously updated
				x_prev = x_curr;
				x_curr = normalize(solve(A_shift, x_curr));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			if (i > max_iter) {
				TH_MATH_ERROR("eigenvalue_rayleigh", i, MathError::NoConvergence);
				return Type(nan());
			}

			return dot(x_curr, transform(A, x_curr));
		}


		/// Compute an eigenvalue of a square matrix and its
		/// corresponding eigenvector (eigenpair) using
		/// the Rayleigh quotient iteration method.
		///
		/// @param A The matrix to compute the eigenvalue of
		/// @param lambda The starting value for the Rayleigh quotient
		/// @param x The starting approximation for the eigenvector
		/// (a random vector is a good choice).
		/// @param v The vector to overwrite with the eigenvector
		/// @param tolerance The minimum difference in norm between subsequent
		/// steps to stop the algorithm at (defaults to ALGEBRA_EIGEN_TOL).
		/// @param max_iter The maximum number of iterations to use
		/// (defaults to ALGEBRA_EIGEN_ITER).
		/// @return An approximate eigenvalue of the matrix, or NaN if the
		/// algorithm did not converge.
		template<typename Matrix, typename Vector1, typename Vector2,
		typename T = matrix_element_t<Matrix>>
		inline auto eigenpair_rayleigh(
			const Matrix& A, const T& lambda, const Vector1& x, Vector2& v,
			real tolerance = ALGEBRA_EIGEN_TOL, unsigned int max_iter = ALGEBRA_EIGEN_ITER) {

			using Type = matrix_element_t<Matrix>;

			if (!is_square(A)) {
				TH_MATH_ERROR("eigenpair_rayleigh", is_square(A), MathError::InvalidArgument);
				return Type(nan());
			}

			if (A.rows() != x.size()) {
				TH_MATH_ERROR("eigenpair_rayleigh", x.size(), MathError::InvalidArgument);
				return Type(nan());
			}

			// Keep track of the shifted matrix
			Matrix A_shift = A;
			mat_shift_diagonal(A_shift, -lambda);

			Type lambda_prev = lambda;
			Type lambda_curr = lambda;

			Vector1 x_prev = normalize(x);
			Vector1 x_curr = normalize(solve(A_shift, x));

			unsigned int i;

			for (i = 1; i <= max_iter; ++i) {
					
				// Update the eigenvalue approximation
				// using Rayleigh quotient (avoiding normalization)
				lambda_prev = lambda_curr;
				lambda_curr = dot(x_curr, transform(A, x_curr));

				// Shift the diagonal by the difference between
				// subsequent eigenvalues steps, to avoid copying matrix A
				mat_shift_diagonal(A_shift, lambda_prev - lambda_curr);

				// Solve the linear system each time, as the shifted
				// matrix is continuously updated
				x_prev = x_curr;
				x_curr = normalize(solve(A_shift, x_curr));

				// Stop the algorithm when |x_k+1 +- x_k| is
				// less then the tolerance in module
				if (norm(x_curr - x_prev) <= tolerance ||
					norm(x_curr + x_prev) <= tolerance)
					break;
			}

			if (i > max_iter) {
				TH_MATH_ERROR("eigenpair_rayleigh", i, MathError::NoConvergence);
				return Type(nan());
			}

			// Overwrite with the eigenvector
			vec_copy(v, x_curr);

			return dot(x_curr, transform(A, x_curr));
		}
	}
}

#endif
