
///
/// @file special.h Special functions
///

#ifndef THEORETICA_SPECIAL_H
#define THEORETICA_SPECIAL_H

#include "constants.h"
#include "real_analysis.h"


namespace theoretica {

	/// @namespace theoretica::special Special functions
	namespace special {


		/// Gamma special function of positive integer argument
		///
		/// @param k The positive integer argument
		/// @result The Gamma function computed using the factorial.
		inline real gamma(unsigned int k) {

			if(k == 0) {
				TH_MATH_ERROR("gamma", k, OUT_OF_DOMAIN);
				return nan();
			}

			return fact(k - 1);
		}


		/// Half Gamma special function, defined as 
		/// HG(n) = Gamma(n / 2)
		/// for any positive integer n.
		///
		/// @param k The positive integer argument
		/// @result The Gamma function computed using
		/// the factorial or double factorial identity.
		inline real half_gamma(unsigned int k) {

			if(k == 0) {
				TH_MATH_ERROR("half_gamma", k, OUT_OF_DOMAIN);
				return nan();
			}

			return (k % 2 == 0)
		        ? fact<uint64_t>(k / 2 - 1)
		        : double_fact<uint64_t>(k - 2) * SQRTPI / (1 << ((k - 1) / 2));
		}


		/// Log Gamma special function of real argument.
		/// This function uses Lanczos' approximation with gamma = 5.
		///
		/// @param x The real argument
		/// @result The logarithm of the Gamma function of x
		inline real lngamma(real x) {

			// Reflection formula for negative values
			if(x < 0) {

				// Check for negative values of Gamma(x)
				if(floor(-x) % 2 == 0) {
					TH_MATH_ERROR("lngamma", x, OUT_OF_DOMAIN);
					return nan();
				}

				return ln(PI / sin(PI * x)) - lngamma(1 - x);
			}

			// Lanczos' coefficients
			const real c[7] = {
				1.000000000178,
				76.180091729400,
				-86.505320327112,
				24.014098222230,
				-1.231739516140,
				0.001208580030,
				-0.000005363820
			};

			// Simplified logarithmic formula
			// for Lanczos' approximation
			real A5 = c[0];

			for (int i = 1; i < 7; ++i)
				A5 += c[i] / (x + i - 1);

			return (x - 0.5) * (ln(x + 4.5) - 1)
					- 5 + ln(SQRTPI * SQRT2 * A5);
		}


		/// Gamma special function of real argument.
		/// This function uses Lanczos' approximation with gamma = 5.
		///
		/// @param x The real argument
		/// @result The Gamma function of x
		inline real gamma(real x) {

			const real x_fract = fract(x);

			// Check if x is a pole or an integer number
			if(x_fract < MACH_EPSILON) {

				if(x <= 0) {
					TH_MATH_ERROR("gamma", x, OUT_OF_DOMAIN);
					return inf();
				} else
					return gamma((unsigned int) x);
			}

			// Check for negative values of Gamma(x)
			// and use the translation identity
			if(x < 0 && floor(-x) % 2 == 0)
				return exp(lngamma(x + 1)) / x;

			// Compute the Gamma function as the exponential
			// of the log Gamma function which uses Lanczos'
			// approximation
			return exp(lngamma(x));
		}


		/// Pi special function of real argument
		///
		/// @param x The real argument
		/// @return The Pi function of x, equal to Gamma(x + 1)
		inline real pi(real x) {
			return gamma(x + 1);
		}


		/// Beta special function of real argument
		///
		/// @param x1 The first real argument
		/// @param x2 The second real argument
		/// @return The Beta function of x1 and x2
		inline real beta(real x1, real x2) {

			return exp(lngamma(x1) + lngamma(x2) - lngamma(x1 + x2));
		}

	}

}

#endif
