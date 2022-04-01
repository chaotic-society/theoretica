
///
/// @file autodiff.h Automatic differentiation
///

#ifndef UROBORO_AUTODIFF_H
#define UROBORO_AUTODIFF_H

#include "dual.h"


namespace uroboro {


	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation
	template<unsigned int N>
	inline vec<N> gradient(dual(*f)(vec<N, dual>), vec<N, real> x) {

		vec<N, real> res;
		vec<N, dual> dual_x;

		for (int i = 0; i < N; ++i)
			dual_x[i] = x[i];

		for (int i = 0; i < N; ++i) {
			dual_x[i].b = 1;
			res.at(i) = f(dual_x).Dual();
			dual_x[i].b = 0;
		}

		return res;
	}

}

#endif
