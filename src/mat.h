#ifndef UROBORO_MATRIX_H
#define UROBORO_MATRIX_H

#include "./constants.h"
#include "./common.h"
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

		inline mat() {
			make_null();
		}

		inline mat(const mat<N, K>& other) {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = other.data[i][l];
				}
			}
		}

		inline mat<N, K>& operator=(const mat<N, K>& other) {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = other.data[i][l];
				}
			}
			return *this;
		}

		inline mat(real diagonal) {
			make_null();
			int diag_n = min(N, K);
			for (int i = 0; i < diag_n; ++i) {
				data[i][i] = diagonal;
			}
		}

		inline ~mat() {}

		inline vec<K> get_column(int l) const {
			vec<K> column;
			for (int i = 0; i < K; ++i) {
				column.data[i] = data[l][i];
			}
			return column;
		}

		inline vec<K> operator[](int l) const {
			return get_column(l);
		}

		inline vec<N> get_row(int l) const {
			vec<N> row;
			for (int i = 0; i < N; ++i) {
				row.data[i] = data[i][l];
			}
			return row;
		}

		inline void set_column(unsigned int l, const vec<K>& column) {
			for (int i = 0; i < K; ++i) {
				data[l][i] = column.data[i];
			}
		}

		inline void set_row(unsigned int l, const vec<N>& row) {
			for (int i = 0; i < N; ++i) {
				data[i][l] = row.data[i];
			}
		}

		inline void make_null() {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = 0;
				}
			}
		}

		inline mat<N, K> operator+(const mat<N, K>& other) const {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] + other.data[i][l];
				}
			}
			return res;
		}

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

		inline vec<K> operator*(const vec<N>& v) const {
			return transform(v);
		}

		template<unsigned int M>
		inline mat<K, M> transform(const mat<M, N>& B) const {

			mat<K, M> res;

			for (int i = 0; i < M; ++i) {
				for (int j = 0; j < K; ++j) {
					res.at(j, i) = get_row(j) * B.get_column(i);
				}
			}

			return res;
		}

		template<unsigned int M>
		inline mat<K, M> operator*(const mat<M, N>& B) const {
			return transform(B);
		}

		// operator += -= *= /= ...

		inline void transpose() {
			if(!is_square())
				// throw
				return;

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

		inline mat<K, N> transposed() const {
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

		inline void set(real a, unsigned int column, unsigned int row) {
			data[column][row] = a;
		}

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
		inline bool is_square() const {
			return N == K;
		}

		inline bool is_diagonal() const {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(i != l && data[i][l] != 0)
						return false;
				}
			}
			return true;
		}

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


		inline static mat<N, K> identity() {
			return mat<N, K>(1.0);
		}

		inline static mat<N, K> diagonal(real diag) {
			return mat<N, K>(diag);
		}


		inline static mat<4, 4> translation(vec3 t) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.at(3, 0) = t[0];
			m.at(3, 1) = t[1];
			m.at(3, 2) = t[2];

			return m;
		}

		inline static mat<4, 4> translation(real x, real y, real z) {

			mat<4, 4> m = mat<4, 4>(1.0);

			m.at(3, 0) = x;
			m.at(3, 1) = y;
			m.at(3, 2) = z;

			return m;
		}


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


		// static mat<N, K> shear(real scale) {
		// 	return diagonal(scale);
		// }

	};

	// Square matrices types
	using mat2 = mat<2, 2>;
	using mat3 = mat<3, 3>;
	using mat4 = mat<4, 4>;

}

#endif
