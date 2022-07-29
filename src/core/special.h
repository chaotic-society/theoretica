
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


		/// Gamma special function of real argument
		///
		/// @param x The real argument
		///
		/// @note This approximation is accurate to 4 significant
		/// digits for positive x, 2 significant digits
		/// for negative x and exact for integer values.
		/// @todo Improve the approximation
		inline real gamma(real x) {

			real x_fract = fract(x);

			// Identity with the factorial
			if(x_fract <= MACH_EPSILON) {

				if(x >= 1) {
					return fact(int(x) - 1);
				} else {
					TH_MATH_ERROR("gamma", x, OUT_OF_DOMAIN);
					return nan();
				}
			}

			real mul = 1;

			// Recursion relation for the Gamma function
			// used for domain reduction to [1, 2]
			while(x > 2) {
				x -= 1;
				mul *= x;
			}

			while(x < 1) {
				mul /= x;
				x += 1;
			}

			// Sixth degree interpolating polynomial in [1, 2]
			return mul * (3.0569 + x * (-4.34693 + x * (3.25067 + x * (-1.12613 + x * 0.165215))));

			// Fourth degree
			// {1.02447174185, 0.9864538837131499882884},
			// {1.20610737385, 0.9165709651955217873784},
			// {1.5, 0.8862269254527580136491},
			// {1.79389262615, 0.9297768595085857482919},
			// {1.97552825815, 0.9898991992254261074602}

			// Eighth degree
			// {1.00759612349, 0.9898991992254261074602},
			// {1.06698729811, 0.9655178178006569003207},
			// {1.17860619516, 0.924134227524792371539},
			// {1.32898992834, 0.8935010083062180359727},
			// {1.5, 0.8862269254527580136491},
			// {1.67101007166, 0.9034652031152693621302},
			// {1.82139380484, 0.9372370730882858255166},
			// {1.93301270189, 0.9735038451293283591075},
			// {1.99240387651, 0.996812206109454852289}
		}


		/// Pi special function of real argument
		///
		/// @param x The real argument
		inline real pi(real x) {
			return gamma(x + 1);
		}


		/// Beta special function of real argument
		///
		/// @param x1 The first real argument
		/// @param x2 The second real argument
		inline real beta(real x1, real x2) {
			return (gamma(x1) * gamma(x2)) / gamma(x1 + x2);
		}

	}

}

#endif
