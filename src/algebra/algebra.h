
///
/// @file algebra.h Linear algebra routines.
/// This file implements all linear algebra routines of the library,
/// working on templated data structures. The Matrix template must 
/// be a class with these methods:
/// - get(i, j)  Get the element in the i-th row and j-th column
/// - at(i, j)  Get a reference to the element in the i-th row and j-th column
/// - size()  Get the number of elements of the matrix or vector (N. rows * N. columns)
/// - rows()  Get the number of rows of the matrix (not defined for vectors)
/// - cols()  Get the number of columns of the matrix (not defined for vectors)
/// - resize() Change or set the size of the matrix
/// The Vector template must be a class with these methods:
/// - get(i)  Get the element in the i-th place
/// - at(i)  Get a reference to the element in the i-th place
/// - size()  Get the total number of elements of the vector
/// - resize()  Change or set the size of the vector
///

#ifndef THEORETICA_ALGEBRA_H
#define THEORETICA_ALGEBRA_H

#include "../complex/complex.h"



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

			using Type = decltype(m.get(0, 0));

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.at(i, j) = (Type) (i == j ? nan() : 0);

			return m;
		}


		/// Overwrite the given vector with the error
		/// vector with NaN values. This function is used
		/// to signal an error.
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_error(Vector& v) {

			using Type = decltype(v.get(0));

			for (unsigned int i = 0; i < v.size(); ++i)
				v.at(i) = (Type) nan();

			return v;
		}


		/// Overwrite a matrix with the identity matrix
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& make_identity(Matrix& m) {

			using Type = decltype(m.get(0, 0));

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.at(i, j) = (Type) (i == j ? 1 : 0);

			return m;
		}


		/// Overwrite a matrix with all zeroes
		/// @param m The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Matrix>
		inline Matrix& mat_zeroes(Matrix& m) {

			using Type = decltype(m.get(0, 0));

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.at(i, j) = (Type) 0;

			return m;
		}


		/// Overwrite a vector with all zeroes
		/// @param v The vector to overwrite
		/// @return A reference to the overwritten vector
		template<typename Vector>
		inline Vector& vec_zeroes(Vector& v) {

			using Type = decltype(v.get(0));

			for (unsigned int i = 0; i < v.size(); ++i)
				v.at(i) = (Type) 0;

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
					dest.at(i, j) = src.get(i, j);

			return dest;
		}


		/// Copy a matrix by overwriting another.
		/// Equivalent to the operation dest = src
		/// @param dest The matrix to overwrite
		/// @param src The matrix to copy
		/// @return A reference to the overwritten matrix
		template<typename Vector1, typename Vector2>
		inline Vector1& vec_copy(Vector1& dest, const Vector2& src) {

			dest.resize(src.size());

			for (unsigned int i = 0; i < src.size(); ++i)
				dest.at(i) = src.get(i);

			return dest;
		}


		/// Returns the square of the Euclidean/Hermitian norm
		/// of the given vector
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto sqr_norm(const Vector& v) {

			using Type = decltype(v.get(0));
			Type sum = (Type) 0;

			// Use conjugation for complex numbers
			if(is_complex_type<Type>::value)
				for (unsigned int i = 0; i < v.size(); ++i)
					sum += v.get(i) * conjugate(v.get(i));
			else
				for (unsigned int i = 0; i < v.size(); ++i)
					sum += v.get(i) * v.get(i);

			return sum;
		}


		/// Returns the Euclidean/Hermitian norm
		/// of the given vector
		/// @param v The vector to compute the norm of
		/// @return The norm of the given vector
		template<typename Vector>
		inline auto norm(const Vector& v) {

			return (decltype(v.get(0))) sqrt(sqr_norm(v));
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
				r.at(i) /= m;

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
				v.at(i) /= m;

			return v;
		}


		/// Computes the dot product between two vectors
		/// @param v1 The first vector
		/// @param v2 The second vector
		/// @return The dot product of the two vectors
		template<typename Vector1, typename Vector2>
		inline auto dot(const Vector1& v1, const Vector2& v2) {

			using Type = decltype(v1.get(0));

			if(v1.size() != v2.size()) {
				TH_MATH_ERROR("algebra::dot", v1.size(), INVALID_ARGUMENT);
				return (Type) nan();
			}

			Type sum = 0;

			// Use conjugation for complex numbers
			if /*constexpr*/ (is_complex_type<Type>::value)
				for (unsigned int i = 0; i < v1.size(); ++i)
					sum += v1.get(i) * conjugate(v2.get(i));
			else
				for (unsigned int i = 0; i < v1.size(); ++i)
					sum += v1.get(i) * v2.get(i);

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

			v3.at(0) = v1.get(1) * v2.get(2) - v1.get(2) * v2.get(1);
			v3.at(1) = v1.get(2) * v2.get(0) - v1.get(0) * v2.get(2);
			v3.at(2) = v1.get(0) * v2.get(1) - v1.get(1) * v2.get(0);

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
					res.at(j, i) = m.get(i, j);

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
					
					const auto buff = m.at(i, j);
					m.at(i, j) = m.at(j, i);
					m.at(j, i) = buff;
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
					dest.at(j, i) = src.get(i, j);

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
					res.at(j, i) = conjugate(m.get(i, j));

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
					
					const auto buff = m.at(i, j);
					m.at(i, j) = conjugate(m.at(j, i));
					m.at(j, i) = conjugate(buff);
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
					dest.at(j, i) = conjugate(src.get(i, j));

			return dest;
		}


		/// Invert the given matrix and overwrite it.
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

			using Type = decltype(src.get(0, 0));

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
				if(A.at(i, i) == (Type) 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < src.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.at(j, i) != (Type) 0) {

							for (unsigned int k = 0; k < src.rows(); ++k) {
								A.at(i, k) += A.at(j, k);
								dest.at(i, k) += dest.at(j, k);
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

				auto inv_pivot = ((Type) 1.0) / A.at(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = 0; j < src.rows(); ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					const auto coeff = A.at(j, i) * inv_pivot;
					
					for (unsigned int k = 0; k < src.rows(); ++k) {
						A.at(j, k) -= coeff * A.at(i, k);
						dest.at(j, k) -= coeff * dest.at(i, k);
					}
				}

				// Divide the current row by the pivot
				for (unsigned int j = 0; j < src.rows(); ++j) {
					A.at(i, j) *= inv_pivot;
					dest.at(i, j) *= inv_pivot;
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

			auto sum = m.get(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				sum += m.get(i, i);

			return sum;
		}


		/// Compute the product of the elements of
		/// the main diagonal of a generic matrix
		/// @param m The input matrix
		/// @return The product of all the elements of the
		/// main diagonal of the input matrix
		template<typename Matrix>
		inline auto diagonal_product(const Matrix& m) {

			auto mul = m.get(0, 0);
			const size_t n = min(m.rows(), m.cols());

			for (unsigned int i = 1; i < n; ++i)
				mul *= m.get(i, i);

			return mul;
		}


		/// Compute the determinant of a square matrix.
		/// Gauss Jordan elimination is used to reduce the
		/// matrix to a triangular matrix.
		/// @param m The matrix to compute the determinant of
		/// @return The determinant of the matrix
		template<typename Matrix>
		inline auto det(const Matrix& m) {

			using Type = decltype(m.get(0, 0));
			
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
				if(A.at(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (unsigned int j = i + 1; j < A.rows(); ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero.
						// The determinant does not change
						// when adding a row to another one
						if(A.at(j, i) != (Type) 0) {

							for (unsigned int k = 0; k < A.rows(); ++k) {
								A.at(i, k) += (Type) A.at(j, k);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						return (Type) 0;
					}
				}

				const auto inv_pivot = ((Type) 1.0) / A.at(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (unsigned int j = i + 1; j < A.rows(); ++j) {

					// Multiplication coefficient for
					// the elision of Ajk
					const auto coeff = (Type) A.at(j, i) * inv_pivot;

					// The coefficient does not change
					// when adding a linear combination
					// of a row to another
					for (unsigned int k = 0; k < A.rows(); ++k) {
						A.at(j, k) -= (Type) A.at(i, k) * coeff;
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
			return m.get(0, 0) * m.get(1, 1) - m.get(1, 0) * m.get(0, 1);
		}


		/// Return the determinant if the matrix is 3x3.
		/// @note No error checking is performed on the matrix size
		template<typename Matrix>
		inline real det_3x3(const Matrix& m) {
			return	m.get(0, 0) * (m.get(1, 1) * m.get(2, 2) - m.get(2, 1) * m.get(1, 2)) -
					m.get(0, 1) * (m.get(1, 0) * m.get(2, 2) - m.get(2, 0) * m.get(1, 2)) +
					m.get(0, 2) * (m.get(1, 0) * m.get(2, 1) - m.get(2, 0) * m.get(1, 1));
		}


		/// Multiply a matrix by a scalar of any compatible type
		/// @param a A scalar value
		/// @param m The matrix to multiply
		/// @return A reference to the multiplied matrix
		template<typename Field, typename Matrix>
		inline Matrix& mat_scalmul(Field a, Matrix& m) {

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					m.at(i, j) *= a;

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
					dest.at(i, j) = a * src.get(i, j);

			return dest;
		}


		/// Multiply a vector by a scalar of any compatible type
		/// @param a A scalar value
		/// @param v The vector to multiply
		/// @return A reference to the multiplied vector
		template<typename Field, typename Vector>
		inline Vector& vec_scalmul(Field a, Vector& v) {

			for (unsigned int i = 0; i < v.size(); ++i)
				v.at(i) *= a;

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
				dest.at(i) = a * src.get(i);

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
					res.at(i) += A.get(i, j) * v.get(j);

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
					res.at(i) += A.get(i, j) * v.get(j);

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
					res.at(i) += A.get(i, j) * v.get(j);

			return res;
		}


		/// Construct a matrix with the given vector as diagonal
		/// and zeroes everywhere else
		/// @param The vector of diagonal elements
		/// @param res The matrix to overwrite
		/// @return A reference to the overwritten matrix
		template<typename Vector, typename Matrix>
		inline Matrix& diagonal(Matrix& res, const Vector& v) {

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
					res.at(i, j) = (i == j) ? v.get(i) : 0;

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
					A.at(i, j) = A.get(i, j) + B.get(i, j);

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
					res.at(i, j) = A.get(i, j) + B.get(i, j);

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
					A.at(i, j) = A.get(i, j) - B.get(i, j);

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
					res.at(i, j) = A.get(i, j) - B.get(i, j);

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
					A.at(i, j) = A.get(i, j) * alpha + B.get(i, j) * beta;

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
					res.at(i, j) = A.get(i, j) * alpha + B.get(i, j) * beta;

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
						res.at(i, j) += A.get(i, k) * B.get(k, j);

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
						res.at(i, j) += A.get(i, k) * B.get(k, j);

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
					if(abs(A.get(i, j) - B.get(i, j)) > tolerance)
						return false;

			return true;
		}


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
					if(i != j && abs(m.get(i, j)) > tolerance)
						return false;

			return true;
		}

		/// Returns whether the matrix is symmetric
		/// @param m The matrix to consider
		/// @param tolerance The tolerance to allow for
		/// in the comparison, defaults to 10 * MACH_EPSILON
		/// @return A boolean value
		template<typename Matrix>
		inline bool is_symmetric(const Matrix& m, real tolerance = 10 * MACH_EPSILON) {

			if(!is_square(m))
				return false;

			for (unsigned int i = 0; i < m.rows(); ++i)
				for (unsigned int j = 0; j < m.cols(); ++j)
					if(abs(m.get(i, j) - m.get(j, i)) > tolerance)
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
				v1.at(i) = v1.get(i) + v2.get(i);

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
				res.at(i) = v1.get(i) + v2.get(i);

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
				v1.at(i) = v1.get(i) - v2.get(i);

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
				res.at(i) = v1.get(i) - v2.get(i);

			return res;
		}


		// Linear transformations


		/// Returns the identity matrix. Size parameters are used
		/// only for dynamically allocated matrix types.
		/// @return The identity matrix of the given type
		template<typename Matrix>
		inline Matrix identity(unsigned int rows = 0, unsigned int cols = 0) {
			
			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);
			
			make_identity(m);
			return m;
		}


		/// Returns the a matrix with the given diagonal. The function
		/// without any parameters is used for statically
		/// allocated matrix types.
		/// @return The diagonal matrix of the given type
		template<typename Matrix, typename Vector>
		inline Matrix diagonal(
			const Vector& v, unsigned int rows = 0, unsigned int cols = 0) {
			
			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);
			
			diagonal(m, v);
			return m;
		}


		/// Returns a translation matrix, that is,
		/// an NxN matrix which describes a translation
		/// in N-1 dimensions.
		/// @param v The translation vector of dimension N-1
		/// @return The translation matrix
		template<typename Matrix, typename Vector>
		inline Matrix translation(
			const Vector& v, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(v.size() != (m.rows() - 1)) {
				TH_MATH_ERROR("algebra::translation", v.size(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			make_identity(m);
				
			// The translation matrix in projective
			// geometry is constructed by setting 
			// the last column's values to the vector,
			// while other elements form the identity
			for (unsigned int i = 0; i < m.rows() - 1; ++i)
				m.at(i, m.cols() - 1) = v[i];

			return m;
		}


		/// Returns a matrix representing a 2D rotation
		/// @param theta The angle of rotation
		/// @return The rotation matrix
		template<typename Matrix>
		inline Matrix rotation_2d(
			real theta, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 2) {
				TH_MATH_ERROR("algebra::rotation_2d", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 2) {
				TH_MATH_ERROR("algebra::rotation_2d", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			const real s = sin(theta);
			const real c = cos(theta);

			if(m.rows() > 2 || m.cols() > 2)
				make_identity(m);

			m.at(0, 0) = c;
			m.at(0, 1) = -s;
			m.at(1, 0) = s;
			m.at(1, 1) = c;

			return m;
		}


		/// Returns a matrix representing a 3D rotation
		/// around a given axis.
		/// @param theta Angle of rotation
		/// @param axis Axis of rotation
		/// @return A matrix representing the rotation
		/// of theta radians around the given axis.
		template<typename Matrix, typename Vector>
		inline Matrix rotation_3d(
			real theta, const Vector& axis, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(axis.size() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d", axis.size(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.rows() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.rows() > 3 || m.cols() > 3)
				make_identity(m);

			const real s = sin(theta);
			const real c = cos(theta);

			const real Rx = (real) axis.get(0);
			const real Ry = (real) axis.get(1);
			const real Rz = (real) axis.get(2);

			const real cm1 = (1 - c);

			m.at(0, 0) = c + Rx * Rx * cm1;
			m.at(0, 1) = Rx * Ry * cm1 - Rz * s;
			m.at(0, 2) = Rx * Rz * cm1 - Ry * s;

			m.at(1, 0) = Ry * Rx * cm1 + Rz * s;
			m.at(1, 1) = c + Ry * Ry * cm1;
			m.at(1, 2) = Ry * Rz * cm1 - Rx * s;

			m.at(2, 0) = Rz * Rx * cm1 - Ry * s;
			m.at(2, 1) = Rz * Ry * cm1 + Rx * s;
			m.at(2, 2) = c + Rz * Rz * cm1;

			return m;
		}


		/// Returns a matrix representing a 3D rotation
		/// around the x axis.
		/// @param theta Angle of rotation
		/// @return A matrix representing the rotation
		/// of theta radians around the x axis.
		template<typename Matrix>
		inline Matrix rotation_3d_xaxis(
			real theta, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_xaxis", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_xaxis", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.rows() > 3 || m.cols() > 3)
				make_identity(m);

			const real s = theoretica::sin(theta);
			const real c = theoretica::cos(theta);

			m.at(0, 0) = 1;
			m.at(1, 1) = c;
			m.at(2, 2) = c;

			m.at(1, 2) = -s;
			m.at(2, 1) = s;

			return m;
		}


		/// Returns a matrix representing a 3D rotation
		/// around the y axis.
		/// @param theta Angle of rotation
		/// @return A matrix representing the rotation
		/// of theta radians around the y axis.
		template<typename Matrix>
		inline Matrix rotation_3d_yaxis(
			real theta, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_yaxis", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_yaxis", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.rows() > 3 || m.cols() > 3)
				make_identity(m);

			const real s = theoretica::sin(theta);
			const real c = theoretica::cos(theta);

			m.at(0, 0) = c;
			m.at(1, 1) = 1;
			m.at(2, 2) = c;
			m.at(3, 3) = 1;

			m.at(0, 2) = s;
			m.at(2, 0) = -s;

			return m;
		}


		/// Returns a matrix representing a 3D rotation
		/// around the z axis.
		/// @param theta Angle of rotation
		/// @return A matrix representing the rotation
		/// of theta radians around the z axis.
		template<typename Matrix>
		inline Matrix rotation_3d_zaxis(
			real theta, unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_zaxis", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 3) {
				TH_MATH_ERROR("algebra::rotation_3d_zaxis", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.rows() > 3 || m.cols() > 3)
				make_identity(m);

			const real s = theoretica::sin(theta);
			const real c = theoretica::cos(theta);

			m.at(0, 0) = c;
			m.at(1, 1) = c;
			m.at(2, 2) = 1;
			m.at(3, 3) = 1;

			m.at(0, 1) = -s;
			m.at(1, 0) = s;

			return m;
		}


		/// Returns a perspective projection matrix
		template<typename Matrix>
		inline Matrix perspective(
			real left, real right, real bottom,
			real top, real near, real far,
			unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 4) {
				TH_MATH_ERROR("algebra::perspective", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 4) {
				TH_MATH_ERROR("algebra::perspective", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			mat_zeroes(m);

			m.at(0, 0)  = 2 * near / (right - left);
			m.at(2, 0)  = (right + left) / (right - left);
			m.at(1, 1)  = 2 * near / (top - bottom);
			m.at(2, 1)  = (top + bottom) / (top - bottom);
			m.at(2, 2) = -(far + near) / (far - near);
			m.at(3, 2) = -(2 * far * near) / (far - near);
			m.at(2, 3) = -1;
			m.at(3, 3) = 0;

			return m;
		}


		/// Returns a perspective projection matrix,
		/// using the Field of View as parameter.
		template<typename Matrix>
		inline Matrix perspective_fov(
			real fov, real aspect, real near, real far,
			unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 4) {
				TH_MATH_ERROR("algebra::perspective_fov", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 4) {
				TH_MATH_ERROR("algebra::perspective_fov", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			const real height = near * tan(radians(fov / 2.f));
			const real width = height * aspect;

			return perspective<Matrix>(-width, width, -height, height, near, far);
		}


		/// Returns an orthogonal projection matrix.
		template<typename Matrix>
		inline Matrix ortho(
			real left, real right, real bottom,
			real top, real near, real far,
			unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(m.rows() < 4) {
				TH_MATH_ERROR("algebra::ortho", m.rows(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			if(m.cols() < 4) {
				TH_MATH_ERROR("algebra::ortho", m.cols(), INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			mat_zeroes(m);

			m.at(0, 0)  = 2 / (right - left);
			m.at(3, 0)  = -(right + left) / (right - left);
			m.at(1, 1)  = 2 / (top - bottom);
			m.at(3, 1)  = -(top + bottom) / (top - bottom);
			m.at(2, 2) = -2 / (far - near);
			m.at(3, 2) = -(far + near) / (far - near);

			return m;
		}


		/// Return a transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		template<
			typename Matrix, typename Vector1,
			typename Vector2, typename Vector3>
		inline Matrix look_at(Vector1 camera, Vector2 target, Vector3 up) {

			// Construct an orthonormal basis

			Vector1 x_axis, y_axis, z_axis;
			x_axis.resize(3);
			y_axis.resize(3);
			z_axis.resize(3);

			// z = target - camera
			vec_diff(z_axis, target, camera);
			make_normalized(z_axis);

			// x = z X up
			vec_copy(x_axis, cross(z_axis, up));
			make_normalized(x_axis);

			// y = x X z
			vec_copy(y_axis, cross(x_axis, z_axis));

			// Negate z_axis to have a right-handed system
			vec_scalmul(-1.0, z_axis);

			// Construct the rotation and translation matrix
			Matrix m;
			m.resize(4, 4);

			m.at(0, 0) = x_axis.get(1);
			m.at(0, 1) = x_axis.get(2);
			m.at(0, 2) = x_axis.get(3);
			m.at(0, 3) = dot(camera, vec_scalmul(-1.0, x_axis));

			m.at(1, 0) = y_axis.get(1);
			m.at(1, 1) = y_axis.get(2);
			m.at(1, 2) = y_axis.get(3);
			m.at(1, 3) = dot(camera, vec_scalmul(-1.0, y_axis));

			m.at(2, 0) = z_axis.get(1);
			m.at(2, 1) = z_axis.get(2);
			m.at(2, 2) = z_axis.get(3);
			m.at(2, 3) = dot(camera, vec_scalmul(-1.0, z_axis));

			m.at(3, 0) = 0;
			m.at(3, 1) = 0;
			m.at(3, 2) = 0;
			m.at(3, 3) = 1;

			return m;
		}


		template<typename Matrix>
		inline Matrix symplectic(unsigned int rows = 0, unsigned int cols = 0) {

			Matrix m;
			if(rows && cols)
				m.resize(rows, cols);

			if(rows != cols || (rows % 2 != 0)) {
				TH_MATH_ERROR("algebra::symplectic", rows, INVALID_ARGUMENT);
				mat_error(m);
				return m;
			}

			const unsigned int half = m.rows() / 2;
			mat_zeroes(m);

			for (unsigned int i = 0; i < half; ++i) {
				m.at(i, i + half) = 1;
				m.at(i + half, i) = -1;
			}

			return m;
		}


		/// Sphere inversion of a point with respect to
		/// a sphere of radius r centered at a point c
		/// @param p The vector to transform
		/// @param c The center of the sphere
		/// @param r The radius of the sphere
		template<typename Vector1, typename Vector2>
		inline Vector1 sphere_inversion(
			const Vector1& p, const Vector2& c = Vector2(0), real r = 1) {

			Vector1 q = p - c;
			return c + q * square(r / q.norm());
		}


		// LU decomposition

		// QR decomposition

		// Eigenvalues

		// Eigenvectors


	}

}

#endif
