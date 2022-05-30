
///
/// @file quasirandom.h Quasi-random sequences
///

#ifndef THEORETICA_QUASIRANDOM_H
#define THEORETICA_QUASIRANDOM_H

#include "../core/real_analysis.h"
#include "../algebra/vec.h"


namespace theoretica {


	/// Weyl quasi-random sequence
	/// @param n The index of the element in the sequence
	/// @param alpha The base of the sequence, defaults to 
	/// the inverse of the Golden Section
	/// 
	/// The Weyl sequence is defined as \f$x_n = {n \alpha}\f$,
	/// where \f${ }\f$ is the fractional part.
	/// @note The alpha argument should be an irrational number.
	real qrand_weyl(unsigned int n, real alpha = INVPHI) {
		return fract(n * alpha);
	}


	/// Weyl quasi-random sequence (computed with recurrence relation)
	/// @param prev The previously computed value
	/// @param alpha The base of the sequence, defaults to 
	/// the inverse of the Golden Section
	///
	///	If no arguments are provided or prev is zero,
	/// the function computes the first element of the
	/// Weyl sequence associated to the parameter alpha.
	/// @see qrand_weyl
	real qrand_weyl_recurr(real prev = 0, real alpha = INVPHI) {

		if(prev == 0) {
			prev = qrand_weyl(1, alpha);
			return prev;
		}

		return fract(prev + alpha);
	}


	/// Weyl quasi-random sequence in N dimensions
	/// @param n The index of the element in the sequence
	/// @param alpha The base of the sequence
	/// 
	/// @note The alpha argument should be an irrational number.
	template<unsigned int N>
	vec<N> qrand_weyl_multi(unsigned int n, real alpha) {
		
		vec<N> r;
		for (unsigned int i = 0; i < N; ++i)
			r[i] = fract(n * pow(alpha, i + 1));

		return r;
	}


	/// Weyl quasi-random sequence in 2 dimensions
	/// @param n The index of the element in the sequence
	/// @param alpha The base of the sequence
	/// 
	/// @note The alpha argument should be an irrational number.
	vec2 qrand_weyl2(unsigned int n, real alpha = 0.7548776662466927) {
		return {fract(n * alpha), fract(n * square(alpha))};
	}


}

#endif
