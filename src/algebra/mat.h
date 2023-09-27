
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
#include "./algebra.h"
#include "../complex/complex.h"

#include <array>


namespace theoretica {

	/// 
	/// @class mat
	/// A generic matrix with real entries
	///
	/// @param N The number of rows
	/// @param K The number of columns
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

		/// Initialize a matrix to all zeroes
		inline mat() {
			algebra::mat_zeroes(*this);
		}


		/// Initialize a matrix from another one
		inline mat(const mat<N, K, Type>& other) {
			algebra::mat_copy(*this, other);
		}


		/// Initialize from an initializer list of the rows
		inline mat(const std::initializer_list<std::array<Type, K>>& rows) {

			if(rows.size() != N) {
				TH_MATH_ERROR("mat::mat(std::array<vec<K>, N>)", rows.size(), INVALID_ARGUMENT);
				*this = mat<N, K, Type>(nan());
				return;
			}

			int i = 0;

			for (const auto& r : rows) {
				for (unsigned int j = 0; j < K; ++j)
					iat(i, j) = r[j];
				i++;
			}
		}


		/// Copy constructor
		inline mat<N, K, Type>& operator=(const mat<N, K, Type>& other) {
			return algebra::mat_copy(*this, other);
		}


		/// Construct a diagonal matrix
		inline mat(Type diagonal) {
			algebra::diagonal(*this, vec<N>(diagonal));
		}


		/// Default destructor
		inline ~mat() = default;


		/// Set all elements to zero
		inline void make_zeroes() {
			algebra::mat_zeroes(*this);
		}


		/// Matrix addition
		inline mat<N, K, Type> operator+(const mat<N, K, Type>& other) const {
			mat<N, K, Type> res;
			return algebra::mat_sum(res, *this, other);
		}


		/// Matrix subtraction
		inline mat<N, K, Type> operator-(const mat<N, K, Type>& other) const {
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

			if(scalar == 0) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(res);
			}
			
			return algebra::mat_scalmul(res, 1.0 / scalar, *this);
		}


		/// Transform a vector v by the matrix
		inline vec<N> transform(const vec<K>& v) const {
			return algebra::transform(*this, v);
		}


		/// Transform a vector v by the matrix
		inline vec<N> operator*(const vec<K>& v) const {
			return transform(v);
		}


		/// Matrix multiplication
		template<unsigned int M>
		inline mat<N, M, Type> transform(const mat<K, M, Type>& B) const {
			mat<N, M, Type> res;
			return algebra::mat_mul(res, *this, B);
		}


		/// Matrix multiplication
		template<unsigned int M>
		inline mat<N, M, Type> operator*(const mat<K, M, Type>& B) const {
			return transform(B);
		}


		/// Matrix addition
		inline mat<N, K, Type>& operator+=(const mat<N, K, Type>& other) {
			return algebra::mat_sum(*this, other);
		}


		/// Matrix subtraction
		inline mat<N, K, Type>& operator-=(const mat<N, K, Type>& other) {
			return algebra::mat_diff(*this, other);
		}


		/// Scalar multiplication
		inline mat<N, K, Type>& operator*=(Type scalar) {
			return algebra::mat_scalmul(scalar, *this);
		}


		/// Scalar division
		inline mat<N, K, Type>& operator/=(Type scalar) {

			if(scalar == 0) {
				TH_MATH_ERROR("mat::operator/", scalar, DIV_BY_ZERO);
				return algebra::mat_error(*this);
			}

			return algebra::mat_scalmul(1.0 / scalar, *this);
		}


		/// Matrix multiplication
		inline mat<N, K, Type>& operator*=(const mat<N, K, Type>& B) {
			return (*this = this->operator*(B));
		}


		/// Transpose the matrix itself
		inline mat<N, K, Type>& transpose() {
			static_assert(
				N == K, "The matrix must be square to be transposed in place.");
			return algebra::make_transposed(*this);
		}


		/// Return the transposed matrix
		inline mat<K, N, Type> transposed() const {
			return algebra::transpose<mat<N, K, Type>, mat<K, N, Type>>(*this);
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
		inline Type& iat(unsigned int i, unsigned int j) {

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
		inline Type& at(unsigned int i, unsigned int j) {

#ifdef THEORETICA_MATRIX_LEXIC
			return iat(j, i);
#else
			return iat(i, j);
#endif
		}


		/// Access element at indices i and j
		inline Type& operator()(unsigned int i, unsigned int j) {
			return at(i, j);
		}


		/// Getters and setters
		inline Type iget(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_ROW_FIRST
			return data[i][j];
#else
			return data[j][i];
#endif
		}


		/// Getters and setters
		inline Type get(unsigned int i, unsigned int j) const {

#ifdef THEORETICA_MATRIX_LEXIC
			return iget(j, i);
#else
			return iget(i, j);
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
		inline bool operator==(const mat<N, K, Type>& other) const {
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
		inline static mat<4, 4> look_at(vec<3> camera, vec<3> target, vec<3> up) {
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

						if(abs(iget(i, j)) < MACH_EPSILON)
							res << "0";
						else
							res << iget(i, j);
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
