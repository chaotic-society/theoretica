
///
/// @file fft.h Fast Fourier Transform
///

#ifndef THEORETICA_FFT_H
#define THEORETICA_FFT_H

#include "../core/bit_op.h"
#include "../algebra/algebra_types.h"
#include "../algebra/algebra.h"
#include "../complex/complex.h"


namespace theoretica {


	/// @namespace theoretica::signal Signal processing module.
	namespace signal {


		/// Compute the Fast Fourier Transform of a set of data points.
		/// Bit reversion is used on the indices to simplify the resulting calculations.
		///
		/// @param x The set of data points in the time domain
		/// @param inverse Whether to run the inverse transform (defaults to false)
		/// @return The data in the frequency domain
		template<typename ReturnVector = cvec, typename InputVector = cvec>
		inline ReturnVector fft(const InputVector& x, bool inverse = false) {

			// Resulting vector in the frequency domain
			ReturnVector k = x;
			const unsigned int N = x.size();

			// Compute the logarithm of the size
			const unsigned int log2N = ilog2(N);

			// Enforce power of 2 vector size
			if (N != (unsigned int) (1 << log2N)) {
				algebra::vec_error(k);
				return k;
			}

			// Bit reverse
			swap_bit_reverse(k, log2N);


			for (unsigned int p = 1; p <= log2N; p++) {

				const unsigned int m = 1 << p;
				const unsigned int offset = m / 2;
				
				complex<real> w (1.0, 0.0);
				const real sign = (inverse ? 1.0 : -1.0);

				// Phase shift between iterations
				const complex<real> phase = complex<real>(
					cos(sign * 2 * PI / m),
					sin(sign * 2 * PI / m)
				);

				for (unsigned int j = 0; j < offset; j++) {

					for (unsigned int i = j; i < N; i += m) {

						const complex<real> t = w * k[i + offset];
						k[i + offset] = k[i] - t;
						k[i] += t;
					}
					
					w *= phase;
				}
			}

			// The normalization constant is 1/N
			if (inverse)
				k /= x.size();

			return k;
		}


		/// Compute the Inverse Fast Fourier Transform of a set of data points.
		/// Bit reversion is used on the indices to simplify the resulting calculations.
		///
		/// @param k The set of data points in the frequency domain
		/// @return The data in the time domain
		template<typename Vector1 = cvec, typename Vector2 = cvec>
		inline Vector2 ifft(const Vector1& k) {
			return fft(k, true);
		}
	}
}


#endif
