
///
/// @file transform.h Linear transformations such as rotations and projective geometry.
///

#ifndef THEORETICA_TRANSFORM_H
#define THEORETICA_TRANSFORM_H

#include "./algebra.h"
#include "../core/error.h"


namespace theoretica {

	namespace algebra {


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


		/// Construct a matrix with the given vector as diagonal
		/// and zeroes everywhere else, overwriting the given matrix.
		///
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
					res(i, j) = (i == j) ? v[i] : 0;

			return res;
		}


		/// Returns the a matrix with the given diagonal. The function
		/// without any parameters is used for statically
		/// allocated matrix types.
		///
		/// @param v The vector of the diagonal elements
		/// @param rows The number of rows of the matrix
		/// @param cols The number of columns of the matrix
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
				m(i, m.cols() - 1) = v[i];

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

			m(0, 0) = c;
			m(0, 1) = -s;
			m(1, 0) = s;
			m(1, 1) = c;

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

			const real Rx = (real) axis[0];
			const real Ry = (real) axis[1];
			const real Rz = (real) axis[2];

			const real cm1 = (1 - c);

			m(0, 0) = c + Rx * Rx * cm1;
			m(0, 1) = Rx * Ry * cm1 - Rz * s;
			m(0, 2) = Rx * Rz * cm1 + Ry * s;

			m(1, 0) = Ry * Rx * cm1 + Rz * s;
			m(1, 1) = c + Ry * Ry * cm1;
			m(1, 2) = Ry * Rz * cm1 - Rx * s;

			m(2, 0) = Rz * Rx * cm1 - Ry * s;
			m(2, 1) = Rz * Ry * cm1 + Rx * s;
			m(2, 2) = c + Rz * Rz * cm1;

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

			const real s = sin(theta);
			const real c = cos(theta);

			m(0, 0) = 1;
			m(1, 1) = c;
			m(2, 2) = c;

			m(1, 2) = -s;
			m(2, 1) = s;

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

			const real s = sin(theta);
			const real c = cos(theta);

			m(0, 0) = c;
			m(1, 1) = 1;
			m(2, 2) = c;
			m(3, 3) = 1;

			m(0, 2) = s;
			m(2, 0) = -s;

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

			const real s = sin(theta);
			const real c = cos(theta);

			m(0, 0) = c;
			m(1, 1) = c;
			m(2, 2) = 1;
			m(3, 3) = 1;

			m(0, 1) = -s;
			m(1, 0) = s;

			return m;
		}


		/// Returns a perspective projection matrix with adjustable view volume boundaries.
		/// @tparam Matrix The type of matrix to be returned.
		/// @param left The left boundary of the view volume.
		/// @param right The right boundary of the view volume.
		/// @param bottom The bottom boundary of the view volume.
		/// @param top The top boundary of the view volume.
		/// @param near The near clipping plane distance.
		/// @param far The far clipping plane distance.
		/// @param rows Optional number of rows to set for the matrix (default is 0).
		/// @param cols Optional number of columns to set for the matrix (default is 0).
		/// @return A perspective projection matrix with dimensions at least 4x4.
		///
		/// This function creates a perspective projection matrix that maps a 3D frustum 
		/// into a 2D plane. If the matrix dimensions are smaller than 4x4, an error is 
		/// triggered. The matrix is initialized to zero, with values set to define the 
		/// specified perspective projection parameters.
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

			m(0, 0)  = 2 * near / (right - left);
			m(2, 0)  = (right + left) / (right - left);
			m(1, 1)  = 2 * near / (top - bottom);
			m(2, 1)  = (top + bottom) / (top - bottom);
			m(2, 2) = -(far + near) / (far - near);
			m(3, 2) = -(2 * far * near) / (far - near);
			m(2, 3) = -1;
			m(3, 3) = 0;

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

			m(0, 0)  = 2 / (right - left);
			m(3, 0)  = -(right + left) / (right - left);
			m(1, 1)  = 2 / (top - bottom);
			m(3, 1)  = -(top + bottom) / (top - bottom);
			m(2, 2) = -2 / (far - near);
			m(3, 2) = -(far + near) / (far - near);

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

			m(0, 0) = x_axis[1];
			m(0, 1) = x_axis[2];
			m(0, 2) = x_axis[3];
			m(0, 3) = dot(camera, vec_scalmul(-1.0, x_axis));

			m(1, 0) = y_axis[1];
			m(1, 1) = y_axis[2];
			m(1, 2) = y_axis[3];
			m(1, 3) = dot(camera, vec_scalmul(-1.0, y_axis));

			m(2, 0) = z_axis[1];
			m(2, 1) = z_axis[2];
			m(2, 2) = z_axis[3];
			m(2, 3) = dot(camera, vec_scalmul(-1.0, z_axis));

			m(3, 0) = 0;
			m(3, 1) = 0;
			m(3, 2) = 0;
			m(3, 3) = 1;

			return m;
		}


		/// Generate a NxN symplectic matrix where N is even.
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
				m(i, i + half) = 1;
				m(i + half, i) = -1;
			}

			return m;
		}


		/// Construct the Hilbert matrix of arbitrary dimension.
		/// The Hilbert matrices are square matrices with particularly high
		/// condition number, which makes them ill-conditioned for numerical calculations.
		/// The elements of the Hilbert matrix are given by
		/// \f$H_{ij} = \frac{1}{i + j - 1}\f$ (for \f$i,j\f$ starting from 1).
		///
		/// @param rows The number of rows (and columns) of the resulting matrix
		/// @return The Hilbert matrix of the given size
		template<typename Matrix>
		inline Matrix hilbert(unsigned int rows = 0) {

			Matrix H;
			if (rows)
				H.resize(rows, rows);

			for (unsigned int i = 0; i < H.rows(); ++i)
				for (unsigned int j = 0; j < H.cols(); ++j)
					H(i, j) = 1.0 / (i + j + 1);

			return H;
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
	}
}

#endif
