#ifndef UROBORO_MATRIX_H
#define UROBORO_MATRIX_H

#include "./constants.h"
#include "./vec.h"

namespace uroboro {

	// N is the number of columns
	// K is the number or rows
	// (column-first order is used for OpenGL)
	template<unsigned int N, unsigned int K>
	class mat {
		public:

		const unsigned int size = N * K;
		const unsigned int column_size = N;
		const unsigned int row_size = K;

		real data[N][K];

		mat() {
			make_null();
		}

		mat(const mat<N, K>& other) {
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					data[i][l] = other.data[i][l];
				}
			}
		}

		mat(real diagonal) {
			make_null();
			int diag_n = min(N, K);
			for (int i = 0; i < diag_n; ++i) {
				data[i][i] = diagonal;
			}
		}

		~mat() {}

		inline vec<K> get_column(int l) {
			vec<K> column;
			for (int i = 0; i < K; ++i) {
				column.data[i] = data[l][i];
			}
			return column;
		}

		inline vec<K> operator[](int l) {
			return get_column(l);
		}

		inline vec<N> get_row(int l) {
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

		inline mat<N, K> operator+(const mat<N, K>& other) {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] + other.data[i][l];
				}
			}
			return res;
		}

		inline mat<N, K> operator-(const mat<N, K>& other) {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] - other.data[i][l];
				}
			}
			return res;
		}

		inline mat<N, K> operator*(real scalar) {
			mat<N, K> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[i][l] = data[i][l] * scalar;
				}
			}
			return res;
		}

		inline vec<K> transform(const vec<N>& v) {
			vec<K> res;
			for (int i = 0; i < K; ++i) {
				res[i] = 0;
				for (int l = 0; l < N; ++l) {
					res[i] += data[i][l] * v.data[l];
				}
			}
			return res;
		}

		inline vec<K> operator*(const vec<N>& v) {
			return transform(v);
		}

		// inline void transpose() {
		// 	mat<N, N> res;
		// 	for (int i = 0; i < N; ++i) {
		// 		for (int l = 0; l < K; ++l) {
		// 			res.data[l][i] = data[i][l];
		// 		}
		// 	}
		// 	mat(res);
		// }

		inline mat<K, N> transposed() {
			mat<K, N> res;
			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					res.data[l][i] = data[i][l];
				}
			}
			return res;
		}

		inline real dot(const vec<N>& v1, const vec<N>& v2) {

			vec<N> o = transform(v2);
			real result = 0;

			for (int i = 0; i < N; ++i) {
				result += v1.data[i] * o.data[i];
			}

			return result;
		}

		inline real& at(unsigned int column, unsigned int row) {
			return data[column][row];
		}

		inline real get(unsigned int column, unsigned int row) {
			return data[column][row];
		}

		inline real set(real a, unsigned int column, unsigned int row) {
			data[column][row] = a;
		}

		inline bool operator==(const mat<N, K>& other) {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(data[i][l] != other.data[i][l])
						return false;
				}
			}

			return true;
		}

		inline bool is_square() {
			return N == K;
		}

		inline bool is_diagonal() {

			for (int i = 0; i < N; ++i) {
				for (int l = 0; l < K; ++l) {
					if(i != l && data[i][l] != 0)
						return false;
				}
			}
			return true;
		}

		inline bool is_symmetric() {

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

	};

}

#endif
