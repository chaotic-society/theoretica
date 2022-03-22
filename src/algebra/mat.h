#ifndef UROBORO_MATRIX_H
#define UROBORO_MATRIX_H

#ifndef UROBORO_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../error.h"
#include "../constants.h"
#include "../real_analysis.h"
#include "./vec.h"


namespace uroboro {

	// KxN matrix implementation
	// N is the number of columns
	// K is the number or rows
	// (column-first order is used for OpenGL)
	template<unsigned int N, unsigned int K>
	class mat {
		public:

		static const unsigned int size = N * K;
		static const unsigned int column_size = N;
		static const unsigned int row_size = K;

		real data[N][K];

		// Initialize a matrix to all zeroes
		inline mat() {
			make_null();
		}

		// Initialize a matrix from another one
		inline mat(const mat<N, K>& other) {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = other.data[i][l];
				}
			}
		}

		// Copy constructor
		inline mat<N, K>& operator=(const mat<N, K>& other) {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = other.data[i][l];
				}
			}
			return *this;
		}

		// Construct a diagonal matrix
		inline mat(real diagonal) {
			make_null();
			int diag_n = min(N, K);
			for (int i = 0; i < diag_n; ++i) {
				data[i][i] = diagonal;
			}
		}

		// Default destructor
		inline ~mat() = default;

		// Get the l-th column of the matrix as a vector
		inline vec<K> get_column(int l) const {
			vec<K> column;
			for (int i = 0; i < K; ++i) {
				column.data[i] = data[l][i];
			}
			return column;
		}

		// Get the l-th column of the matrix as a vector
		inline vec<K> operator[](int l) const {
			return get_column(l);
		}

		// Get the l-th row of the matrix as a vector
		inline vec<N> get_row(int l) const {
			vec<N> row;
			for (int i = 0; i < N; ++i) {
				row.data[i] = data[i][l];
			}
			return row;
		}

		// Se the l-th column of the matrix from a vector
		inline void set_column(unsigned int l, const vec<K>& column) {
			for (int i = 0; i < K; ++i) {
				data[l][i] = column.data[i];
			}
		}

		// Se the l-th row of the matrix from a vector
		inline void set_row(unsigned int l, const vec<N>& row) {
			for (int i = 0; i < N; ++i) {
				data[i][l] = row.data[i];
			}
		}

		// Set all elements to zero
		inline void make_null() {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = 0;
				}
			}
		}

		// Matrix addition
		inline mat<N, K> operator+(const mat<N, K>& other) const {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] + other.data[i][l];
				}
			}
			return res;
		}

		// Matrix subtraction
		inline mat<N, K> operator-(const mat<N, K>& other) const {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] - other.data[i][l];
				}
			}
			return res;
		}

		// Scalar multiplication
		inline mat<N, K> operator*(real scalar) const {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] * scalar;
				}
			}
			return res;
		}

		// Scalar division
		inline mat<N, K> operator/(real scalar) const {

			if(scalar == 0) {
				UMATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return mat<N, K>(nan());
			}

			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] / scalar;
				}
			}
			return res;
		}

		// Transform a vector v by the matrix
		inline vec<K> transform(const vec<N>& v) const {
			vec<K> res;
			for (int i = 0; i < K; ++i) {
				res[i] = 0;
				for (int l = 0; l < N; ++l) {
					res[i] += data[i][l] * v.data[l];
				}
			}
			return res;
		}

		// Transform a vector v by the matrix
		inline vec<K> operator*(const vec<N>& v) const {
			return transform(v);
		}

		// Matrix multiplication
		template<unsigned int M>
		inline mat<K, M> transform(const mat<M, N>& B) const {

			mat<K, M> res;

			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < K; ++j) {
					for (int k = 0; k < N; ++k) {
						res.at(j, i) += data[k][j] * B.data[i][k];
					}
				}
			}

			return res;
		}

		// Matrix multiplication
		template<unsigned int M>
		inline mat<K, M> operator*(const mat<M, N>& B) const {
			return transform(B);
		}


		// Matrix addition
		inline mat<N, K>& operator+=(const mat<N, K>& other) {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] += other.data[i][l];
				}
			}

			return *this;
		}

		// Matrix subtraction
		inline mat<N, K>& operator-=(const mat<N, K>& other) {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] -= other.data[i][l];
				}
			}

			return *this;
		}

		// Scalar multiplication
		inline mat<N, K>& operator*=(real scalar) {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] *= scalar;
				}
			}

			return *this;
		}

		// Scalar division
		inline mat<N, K>& operator/=(real scalar) {

			if(scalar == 0) {
				UMATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return mat<N, K>(nan());
			}

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] /= scalar;
				}
			}

			return *this;
		}

		// Matrix multiplication
		inline mat<N, K>& operator*=(const mat<N, K>& B) {
			return (*this = this->operator*(B));
		}


		// Matrix transposition
		inline void transpose() {

			if(!is_square()) {
				UMATH_ERROR("mat::transpose", K, IMPOSSIBLE_OPERATION);
				// Set all elements to nan ?
				return;
			}

			mat<K, N> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[l][i] = data[i][l];
				}
			}

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = res.data[i][l];
				}
			}
		}

		// Return the transposed matrix
		inline mat<K, N> transposed() const {

			if(!is_square()) {
				UMATH_ERROR("mat::transposed", K, IMPOSSIBLE_OPERATION);
				return diagonal(nan());
			}

			mat<K, N> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[l][i] = data[i][l];
				}
			}
			return res;
		}


		// Calculate the dot product between v1 and v2
		// using this matrix as the product matrix
		inline real dot(const vec<N>& v1, const vec<N>& v2) const {

			vec<N> o = transform(v2);
			real result = 0;

			for (int i = 0; i < N; ++i) {
				result += v1.data[i] * o.data[i];
			}

			return result;
		}


		// Access element at <column, row>
		inline real& at(unsigned int column, unsigned int row) {
			return data[column][row];
		}

		// Getters and setters
		inline real get(unsigned int column, unsigned int row) const {
			return data[column][row];
		}

		inline void set(unsigned int column, unsigned int row, real a) {
			data[column][row] = a;
		}

		// Check whether two matrices are equal element by element
		inline bool operator==(const mat<N, K>& other) const {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(data[i][l] != other.data[i][l])
						return false;
				}
			}

			return true;
		}

		// Matrix types

		// Return whether the matrix is square
		inline bool is_square() const {
			return N == K;
		}

		// Return whether the matrix is diagonal
		inline bool is_diagonal() const {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(i != l && data[i][l] != 0)
						return false;
				}
			}
			return true;
		}

		// Return whether the matrix is symmetric
		inline bool is_symmetric() const {

			if(!is_square())
				return false;

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(i != l && data[i][l] != data[l][i])
						return false;
				}
			}
			return true;
		}


		// Compute the trace (sum of elements on the diagonal) of a square matrix
		inline real trace() {

			if(!is_square()) {
				UMATH_ERROR("mat::trace", K, IMPOSSIBLE_OPERATION);
				return nan();
			}

			real res = 0;

			for (int i = 0; i < N; ++i)
				res += data[i][i];

			return res;
		}


		// Compute the product of the diagonal elements of a square matrix
		inline real diagonal_product() {

			real res = data[0][0];

			for (int i = 1; i < min(N, K); ++i)
				res *= data[i][i];

			return res;	
		}


		// Return the determinant if the matrix is 2x2 (NO error checking)
		inline real det_2x2() {
			return data[0][0] * data[1][1] - data[0][1] * data[1][0];
		}


		// Return the determinant if the matrix is 3x3 (NO error checking)
		inline real det_3x3() {
			return	data[0][0] * (data[1][1] * data[2][2] - data[1][2] * data[2][1]) +
					data[1][0] * (data[0][1] * data[2][2] - data[1][2] * data[2][1]) +
					data[2][0] * (data[0][1] * data[1][2] - data[0][2] * data[1][1]);
		}


		// Compute the determinant of the matrix using
		// Gauss-Jordan lower triangularization
		inline real det_gj() {

			if(!is_square()) {
				UMATH_ERROR("mat::det_gj", N, IMPOSSIBLE_OPERATION);
				return nan();
			}

			mat<N, K> A = mat<N, K>(*this);

			// Iterate on all columns
			for (int i = 0; i < N; ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.at(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (int j = i + 1; j < N; ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero.
						// The determinant does not change
						// when adding a row to another
						if(A.at(i, j) != 0) {

							for (int k = 0; k < N; ++k) {
								A.at(k, i) += A.at(k, j);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						return 0;
					}
				}

				real inv_pivot = 1.0 / A.at(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (int j = i + 1; j < N; ++j) {

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.at(i, j) * inv_pivot;

					// The coefficient does not change
					// when adding a linear combination
					// of a row to another
					for (int k = 0; k < N; ++k) {
						A.at(k, j) -= coeff * A.at(k, i);
					}
				}
			}

			// The determinant of a (lower) triangular matrix
			// is the product of the elements on its diagonal
			return A.diagonal_product();
		}


		// Compute the determinant of the matrix
		inline real det() {

			if(!is_square()) {
				UMATH_ERROR("mat::det", K, IMPOSSIBLE_OPERATION);
				return nan();
			}
			
			if(N == 2)
				return det_2x2();

			if(N == 3)
				return det_3x3();

			return det_gj();
		}


		// Matrix inversion

		// Compute the inverse of a 2x2 matrix (NO ERROR CHECKING)
		inline mat<N, K> inverse_2x2() {

			mat<N, K> B;

			// Exchange elements on the diagonal
			B.at(0, 0) = get(1, 1);
			B.at(1, 1) = get(0, 0);

			// Change sign of the other elements
			B.at(0, 1) = -get(0, 1);
			B.at(1, 0) = -get(1, 0);

			return B / B.det();
		}


		// Compute the inverse of a generic square matrix
		inline mat<N, K> inverse() {

			if(!is_square()) {
				UMATH_ERROR("mat::inverse", N, IMPOSSIBLE_OPERATION);
				return mat<N, K>(nan());
			}

			if(N == 2) {
				return inverse_2x2();
			}

			// Initialize the needed matrices
			// (A|B) is the augmented matrix
			mat<N, K> A = mat<N, K>(*this);
			mat<N, K> B = mat<N, K>::identity();

			// Gauss-Jordan elimination

			// Iterate on all columns
			for (int i = 0; i < N; ++i) {
				
				// Make sure the element on the diagonal
				// is non-zero by adding the first non-zero row
				if(A.at(i, i) == 0) {

					bool flag = false;

					// Iterate on all rows
					for (int j = i + 1; j < N; ++j) {

						// Add the j-th row to the i-th row
						// if Aji is non-zero
						if(A.at(i, j) != 0) {

							for (int k = 0; k < N; ++k) {
								A.at(k, i) += A.at(k, j);
								B.at(k, i) += B.at(k, j);
							}

							flag = true;
							break;
						}
					}

					if(!flag) {
						UMATH_ERROR("mat::inverse", flag, IMPOSSIBLE_OPERATION);
						return mat<N, K>(nan());
					}
				}

				real inv_pivot = 1.0 / A.at(i, i);

				// Use the current row to make all other
				// elements of the column equal to zero
				for (int j = 0; j < N; ++j) {

					// Skip the current row
					if(j == i)
						continue;

					// Multiplication coefficient for
					// the elision of Ajk
					real coeff = A.at(i, j) * inv_pivot;
					
					for (int k = 0; k < N; ++k) {
						A.at(k, j) -= coeff * A.at(k, i);
						B.at(k, j) -= coeff * B.at(k, i);
					}
				}

				// Divide the current row by the pivot
				for (int j = 0; j < N; ++j) {
					A.at(j, i) *= inv_pivot;
					B.at(j, i) *= inv_pivot;
				}
				
			}

			return B;
		}



		// Transformation matrices

		// Get the identity matrix
		inline static mat<N, K> identity() {
			return mat<N, K>(1.0);
		}

		// Get a diagonal matrix
		inline static mat<N, K> diagonal(real diag) {
			return mat<N, K>(diag);
		}

		// Get a 4x4 matrix which translates by {x, y, z}
		inline static mat<4, 4> translation(vec3 t) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.at(3, 0) = t[0];
			m.at(3, 1) = t[1];
			m.at(3, 2) = t[2];

			return m;
		}

		// Get a 4x4 matrix which translates by {x, y, z}
		inline static mat<4, 4> translation(real x, real y, real z) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.at(3, 0) = x;
			m.at(3, 1) = y;
			m.at(3, 2) = z;

			return m;
		}

		// Get a matrix which rotates the 2D plane of <theta> radians
		inline static mat<2, 2> rotation_2x2(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<2, 2> res;
			res.at(0, 0) = c;
			res.at(1, 0) = -s;
			res.at(0, 1) = s;
			res.at(1, 1) = c;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<4, 4> rotation_x_4x4(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.at(0, 0) = 1;
			res.at(1, 1) = c;
			res.at(2, 2) = c;
			res.at(3, 3) = 1;

			res.at(2, 1) = -s;
			res.at(1, 2) = s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<3, 3> rotation_x_3x3(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.at(0, 0) = 1;
			res.at(1, 1) = c;
			res.at(2, 2) = c;

			res.at(2, 1) = -s;
			res.at(1, 2) = s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<4, 4> rotation_y_4x4(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.at(0, 0) = c;
			res.at(1, 1) = 1;
			res.at(2, 2) = c;
			res.at(3, 3) = 1;

			res.at(2, 0) = s;
			res.at(0, 2) = -s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<3, 3> rotation_y_3x3(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.at(0, 0) = c;
			res.at(1, 1) = 1;
			res.at(2, 2) = c;

			res.at(2, 0) = s;
			res.at(0, 2) = -s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<4, 4> rotation_z_4x4(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<4, 4> res = mat<4, 4>(1.0);
			res.at(0, 0) = c;
			res.at(1, 1) = c;
			res.at(2, 2) = 1;
			res.at(3, 3) = 1;

			res.at(1, 0) = -s;
			res.at(0, 1) = s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<3, 3> rotation_z_3x3(real theta) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			mat<3, 3> res = mat<3, 3>(1.0);
			res.at(0, 0) = c;
			res.at(1, 1) = c;
			res.at(2, 2) = 1;

			res.at(1, 0) = -s;
			res.at(0, 1) = s;

			return res;
		}

		// Get a matrix which rotates <theta> radians around the <axis> axis
		inline static mat<4, 4> rotation_4x4(real theta, const vec3& axis) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			real Rx = axis.get(0);
			real Ry = axis.get(1);
			real Rz = axis.get(2);

			real cm1 = (1 - c);

			mat<4, 4> m;

			m.at(0, 0) = c + square(Rx) * cm1;
			m.at(1, 0) = Rx * Ry * cm1 - Rz * s;
			m.at(2, 0) = Rx * Rz * cm1 - Ry * s;
			m.at(3, 0) = 0;

			m.at(0, 1) = Ry * Rx * cm1 + Rz * s;
			m.at(1, 1) = c + square(Ry) * cm1;
			m.at(2, 1) = Ry * Rz * cm1 - Rx * s;
			m.at(3, 1) = 0;

			m.at(0, 2) = Rz * Rx * cm1 - Ry * s;
			m.at(1, 2) = Rz * Ry * cm1 + Rx * s;
			m.at(2, 2) = c + square(Rz) * cm1;
			m.at(3, 2) = 0;

			m.at(0, 3) = 0;
			m.at(1, 3) = 0;
			m.at(2, 3) = 0;
			m.at(3, 3) = 1;

			return m;
		}


		// Get a matrix which rotates <theta> radians around the <axis> axis
		inline static mat<3, 3> rotation_3x3(real theta, const vec3& axis) {

			real s = uroboro::sin(theta);
			real c = uroboro::cos(theta);

			real Rx = axis.get(0);
			real Ry = axis.get(1);
			real Rz = axis.get(2);

			real cm1 = (1 - c);

			mat<3, 3> m;

			m.at(0, 0) = c + square(Rx) * cm1;
			m.at(1, 0) = Rx * Ry * cm1 - Rz * s;
			m.at(2, 0) = Rx * Rz * cm1 - Ry * s;

			m.at(0, 1) = Ry * Rx * cm1 + Rz * s;
			m.at(1, 1) = c + square(Ry) * cm1;
			m.at(2, 1) = Ry * Rz * cm1 - Rx * s;

			m.at(0, 2) = Rz * Rx * cm1 - Ry * s;
			m.at(1, 2) = Rz * Ry * cm1 + Rx * s;
			m.at(2, 2) = c + square(Rz) * cm1;

			return m;
		}


		// Get a scaling matrix by factors <x>, <y> and <z>
		inline static mat<4, 4> scaling_4x4(real x, real y, real z) {

			mat<4, 4> res;
			res.make_null();

			res.at(0, 0) = x;
			res.at(1, 1) = y;
			res.at(2, 2) = z;
			res.at(3, 3) = 1;

			return res;
		}


		// Get a scaling matrix by factors <x>, <y> and <z>
		inline static mat<3, 3> scaling_3x3(real x, real y, real z) {

			mat<3, 3> res;
			res.make_null();

			res.at(0, 0) = x;
			res.at(1, 1) = y;
			res.at(2, 2) = z;

			return res;
		}

		// Get a scaling matrix by <v> factors
		template<unsigned int M>
		inline static mat<N, K> scaling(vec<M> v) {

			mat<N, K> res = mat<N, K>(1.0);

			for (int i = 0; i < min(min(M, N), K); ++i)
				res.at(i, i) = v.get(i);

			return res;
		}


		inline static mat<4, 4> perspective(real left, real right, real bottom, real top, real near, real far) {

			mat<4, 4> result;
			result.make_null();

			result.at(0, 0)  = 2 * near / (right - left);
			result.at(0, 2)  = (right + left) / (right - left);
			result.at(1, 1)  = 2 * near / (top - bottom);
			result.at(1, 2)  = (top + bottom) / (top - bottom);
			result.at(2, 2) = -(far + near) / (far - near);
			result.at(2, 3) = -(2 * far * near) / (far - near);
			result.at(3, 2) = -1;
			result.at(3, 3) = 0;

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

			result.at(0, 0)  = 2 / (right - left);
			result.at(0, 3)  = -(right + left) / (right - left);
			result.at(1, 1)  = 2 / (top - bottom);
			result.at(1, 3)  = -(top + bottom) / (top - bottom);
			result.at(2, 2) = -2 / (far - near);
			result.at(2, 3) = -(far + near) / (far - near);

			return result;
		}

		// Return a transformation matrix that points the
		// field of view towards a given point from the <camera> point
		inline static mat<4, 4> lookAt(vec3 camera, vec3 target, vec3 up) {

			// Construct an orthonormal basis
			vec3 z_axis = (target - camera).normalized();
			vec3 x_axis = cross(z_axis, up).normalized();
			vec3 y_axis = cross(x_axis, z_axis); // Already normalized

			// Negate z_axis to have a right-handed system
			z_axis = -z_axis;

			// Construct the rotation and translation matrix
			mat<4, 4> res;

			res.at(0, 0) = x_axis[0];
			res.at(1, 0) = x_axis[1];
			res.at(2, 0) = x_axis[2];
			res.at(3, 0) = -dot(x_axis, camera);

			res.at(0, 1) = y_axis[0];
			res.at(1, 1) = y_axis[1];
			res.at(2, 1) = y_axis[2];
			res.at(3, 1) = -dot(y_axis, camera);

			res.at(0, 2) = z_axis[0];
			res.at(1, 2) = z_axis[1];
			res.at(2, 2) = z_axis[2];
			res.at(3, 2) = -dot(z_axis, camera);

			res.at(0, 3) = 0;
			res.at(1, 3) = 0;
			res.at(2, 3) = 0;
			res.at(3, 3) = 1;

			return res;
		}

#ifndef UROBORO_NO_PRINT

			// Convert the matrix to string representation
			inline std::string to_string(std::string separator = ", ") const {

				std::stringstream res;

				for (int i = 0; i < row_size; ++i)
					res << get_row(i).to_string(separator) << std::endl;

				return res.str();
			}


			// Stream the matrix in string representation
			// to an output stream (std::ostream)
			friend std::ostream& operator<<(std::ostream& out, const mat<N, K>& obj) {
				return out << obj.to_string();
			}

#endif

	};

	// Square matrices types
	using mat2 = mat<2, 2>;
	using mat3 = mat<3, 3>;
	using mat4 = mat<4, 4>;

}

#endif
