
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

#include "../core/error.h"
#include "../core/constants.h"
#include "../core/real_analysis.h"
#include "./algebra.h"
#include "./vec.h"
#include "../complex/complex.h"


namespace theoretica {

	/// 
	/// @class mat
	/// A generic matrix with real entries
	///
	/// @param N The number of rows
	/// @param K The number of columns
	/// @param Type The type of the elements
	///
	template<unsigned int N, unsigned int K, typename Type = real>
	class mat {
		public:

		static_assert(N > 0, "N cannot be zero.");
		static_assert(K > 0, "K cannot be zero.");

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
		inline mat(const std::initializer_list<std::array<Type, K>>& rows) {

			if(rows.size() != N) {
				TH_MATH_ERROR("mat::mat", rows.size(), INVALID_ARGUMENT);
				algebra::mat_error(*this);
				return;
			}

			int i = 0;
			for (const auto& r : rows) {
				for (unsigned int j = 0; j < K; ++j)
					at(i, j) = r[j];
				i++;
			}
		}


		/// Copy constructor
		template<typename Matrix>
		inline mat<N, K, Type>& operator=(const Matrix& other) {
			return algebra::mat_copy(*this, other);
		}


		/// Construct a diagonal matrix with all equal entries
		mat(Type diagonal) {
			algebra::mat_zeroes(*this);
			for (unsigned int i = 0; i < min(N, K); ++i)
				data[i][i] = diagonal;
		}


		/// Set all elements to zero
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Get the null matrix
		inline static mat<N, K, Type> zeroes() {
			mat<N, K, Type> res;
			algebra::mat_zeroes(res);
			return res;
		}


		/// Matrix addition
		template<typename Matrix>
		inline mat<N, K, Type> operator+(const Matrix& other) const {
			mat<N, K, Type> res;
			return algebra::mat_sum(res, *this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<N, K, Type> operator-(const Matrix& other) const {
			mat<N, K, Type> res;
			return algebra::mat_diff(res, *this, other);
		}


		/// Scalar multiplication
		inline mat<N, K, Type> operator*(Type scalar) const {
			mat<N, K, Type> res;
			return algebra::mat_scalmul(res, scalar, *this);
		}


		/// Friend operator to enable equations of the form
		/// (T) * (mat)
		inline friend mat<N, K, Type> operator*(Type a, const mat<N, K, Type>& B) {
			return B * a;
		}


		/// Scalar division
		inline mat<N, K, Type> operator/(Type scalar) const {

			mat<N, K, Type> res;

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transform a vector v by the matrix
		template<typename Vector>
		inline Vector transform(const Vector& v) const {

			if(v.size() != N) {
				TH_MATH_ERROR("mat::transform", v.size(), INVALID_ARGUMENT);
				Vector res;
				res.resize(N);
				algebra::vec_error(res);
				return res;
			}

			return algebra::transform(*this, v);
		}


		/// Transform a vector by the matrix
		inline vec<N, Type> transform(const vec<K, Type>& v) const {
			return algebra::transform(*this, v);
		}


		/// Transform a vector by the matrix
		inline vec<N, Type> operator*(const vec<K, Type>& v) const {
			return transform(v);
		}


		/// Matrix multiplication
		template<unsigned int M>
		inline mat<N, M, Type> mul(const mat<K, M, Type>& B) const {
			mat<N, M, Type> res;
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
		inline mat<N, K, Type>& operator+=(const Matrix& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction
		template<typename Matrix>
		inline mat<N, K, Type>& operator-=(const Matrix& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication
		inline mat<N, K, Type>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division
		inline mat<N, K, Type>& operator/=(Type scalar) {

			if(abs(scalar) < MACH_EPSILON) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}


		/// Matrix multiplication
		template<typename Matrix>
		inline mat<N, K, Type>& operator*=(const Matrix& B) {
			return (*this = this->operator*(B));
		}


		/// Transpose the matrix itself
		inline mat<N, K, Type>& transpose() {
			static_assert(
				N == K, "The matrix must be square to be transposed in place.");
			return algebra::make_transposed(*this);
		}


		/// Return the transposed matrix, without modifying the
		/// matrix itself.
		inline mat<K, N, Type> transposed() const {
			return algebra::transpose<mat<N, K, Type>, mat<K, N, Type>>(*this);
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
		inline Type& operator()(unsigned int i, unsigned int j) {
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
		inline mat<N, K, Type> inverse() const {
			static_assert(N == K, "The matrix must be square to be invertible.");
			return algebra::inverse(*this);
		}


		/// Invert a generic square matrix
		inline mat<N, K, Type>& invert() {
			static_assert(N == K, "The matrix must be square to be invertible.");
			return algebra::invert(*this);
		}


		/// Compatibility function to allow for allocation
		/// or resizing of dynamic matrices. Since statically
		/// allocated matrices cannot change size, this function
		/// only checks whether the target size is the same
		/// as the matrix's.
		inline void resize(unsigned int n, unsigned int k) const {

			if(N != n) {
				TH_MATH_ERROR("mat::resize", n, INVALID_ARGUMENT);
			} else if(K != k) {
				TH_MATH_ERROR("mat::resize", k, INVALID_ARGUMENT);
			}
		}


		// Transformation matrices


		/// Get the identity matrix
		inline static mat<N, K, Type> identity() {
			return algebra::identity<mat<N, K, Type>>();
		}


		/// Get a diagonal matrix
		inline static mat<N, K, Type> diagonal(Type diag) {
			return mat<N, K, Type>(diag);
		}


		/// Get a 4x4 matrix which translates by {x, y, z}
		template<typename Vector = vec<N - 1>>
		inline static mat<N, K, Type> translation(Vector&& t) {
			return algebra::translation<mat<N, K, Type>>(t);
		}


		/// Get a matrix which rotates the 2D plane of <theta> radians
		inline static mat<N, K, Type> rotation_2d(real theta) {
			static_assert(N >= 2 && K >= 2, "The matrix must be 2x2 or bigger");
			return algebra::rotation_2d<mat<N, K, Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the x axis
		inline static mat<N, K, Type> rotation_3d_xaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_xaxis<mat<N, K, Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the y axis
		inline static mat<N, K, Type> rotation_3d_yaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_yaxis<mat<N, K, Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the z axis
		inline static mat<N, K, Type> rotation_3d_zaxis(real theta) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d_zaxis<mat<N, K, Type>>(theta);
		}


		/// Get a matrix which rotates <theta> radians around the <axis> axis
		template<typename Vector = vec<3>>
		inline static mat<N, K, Type> rotation_3d(real theta, Vector&& axis) {
			static_assert(N >= 3 && K >= 3, "The matrix must be 3x3 or bigger");
			return algebra::rotation_3d<mat<N, K, Type>>(theta, axis);
		}


		inline static mat<N, K, Type> perspective(
			real left, real right, real bottom,
			real top, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective<mat<N, K, Type>>(
				left, right, bottom, top, near, far);
		}


		inline static mat<N, K, Type> perspective_fov(
			real fov, real aspect, real near, real far) {

			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::perspective_fov<mat<N, K, Type>>(fov, aspect, near, far);
		}


		inline static mat<N, K, Type> ortho(
			real left, real right, real bottom, real top, real near, real far) {
			static_assert(N >= 4 && K >= 4, "The matrix must be 4x4 or bigger");
			return algebra::ortho<mat<N, K, Type>>(left, right, bottom, top, near, far);
		}


		/// Return a 4x4 transformation matrix that points the
		/// field of view towards a given point from the <camera> point
		template<typename Vector1, typename Vector2, typename Vector3>
		inline static mat<4, 4> look_at(
			const Vector1& camera, const Vector2& target, const Vector3& up) {
			return algebra::look_at<mat<4, 4>>(camera, target, up);
		}


		/// A symplectic NxN matrix, where \f$N = 2K\f$ for some natural K
		inline static mat<N, K, Type> symplectic() {
			static_assert(N == K && (N % 2 == 0),
				"N must equal K and they should be a multiple of 2");
			return algebra::symplectic<mat<N, K, Type>>();
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


			/// Stream the matrix in string representation to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const mat<N, K, Type>& obj) {
				return out << obj.to_string();
			}

#endif

	};

}

#endif
