
///
/// @file autodiff.h Differential operators using automatic differentiation
///

#ifndef UROBORO_AUTODIFF_H
#define UROBORO_AUTODIFF_H

#include "dual.h"
#include "multidual.h"
#include "dual2.h"

namespace uroboro {


	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	template<unsigned int N>
	inline vec<N> gradient(multidual<N>(*f)(vec<N, multidual<N>>), const vec<N, real>& x) {
		return f(multidual<N>::pack_function_arg(x)).Dual();
	}


	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	///
	/// @note The multidual implementation is more efficient as it does not
	/// need to compute the function value N times and should be preferred.
	/// 
	/// The `mono` suffix is used to emphasize the difference between simple
	/// dual numbers and multidual numbers and to avoid differentiation
	/// between overloads on the user's side.
	template<unsigned int N>
	inline vec<N> gradient_mono(dual(*f)(vec<N, dual>), const vec<N, real>& x) {

		vec<N, real> res;
		vec<N, dual> dual_x;

		for (int i = 0; i < N; ++i)
			dual_x[i] = x.get(i);

		for (int i = 0; i < N; ++i) {
			dual_x[i].b = 1;
			res.at(i) = f(dual_x).Dual();
			dual_x[i].b = 0;
		}

		return res;
	}


	/// Compute the divergence for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation
	template<unsigned int N>
	inline real divergence(multidual<N>(*f)(vec<N, multidual<N>>), const vec<N, real>& x) {

		multidual<N> d = f(multidual<N>::pack_function_arg(x));

		real div = 0;
		for (int i = 0; i < N; ++i)
			div += d.v.at(i);

		return div;
	}


	/// Compute the divergence for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	///
	/// @note The multidual implementation is more efficient as it does not
	/// need to compute the function value N times and should be preferred.
	/// 
	/// The `mono` suffix is used to emphasize the difference between simple
	/// dual numbers and multidual numbers and to avoid differentiation
	/// between overloads on the user's side.
	template<unsigned int N>
	inline real divergence_mono(dual(*f)(vec<N, dual>), const vec<N, real>& x) {

		real res = 0;
		vec<N, dual> dual_x;

		for (int i = 0; i < N; ++i)
			dual_x[i] = x.get(i);

		for (int i = 0; i < N; ++i) {
			dual_x[i].b = 1;
			res += f(dual_x).Dual();
			dual_x[i].b = 0;
		}

		return res;
	}


	/// Compute the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$
	template<unsigned int N, unsigned int M>
	inline mat<M, N> jacobian(vec<M, multidual<N>>(*f)(vec<N, multidual<N>>), vec<N, real> x) {

		vec<M, multidual<N>> res = f(multidual<N>::pack_function_arg(x));

		// Construct the jacobian matrix
		mat<M, N> J;
		for (int i = 0; i < N; ++i) {
			for (int j = 0; j < M; ++j) {
				J.iat(i, j) = res.at(j).Dual().at(i);
			}
		}

		return J;
	}


	/// Compute the directional derivative of a generic function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// @param f The function to partially differentiate
	/// @param x The point to compute the derivative at
	/// @param v The direction to compute the derivative on
	///
	/// \note In most applications, the vector v should be a unit vector,
	/// but the function does not control whether the vector has
	/// unit length or not.
	template<unsigned int N>
	inline vec<N, real> directional_derivative(multidual<N>(*f)(vec<N, multidual<N>>),
		const vec<N, real>& x, const vec<N, real>& v) {
		return v.normalized() * dot(v, gradient(f, x));
	}


	/// Compute the laplacian differential operator for a generic
	/// function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	template<unsigned int N>
	inline real laplacian(dual2(*f)(vec<N, dual2>), const vec<N>& x) {

		real res = 0;
		vec<N, dual2> dual_x;

		for (int i = 0; i < N; ++i)
			dual_x[i] = x.get(i);

		for (int i = 0; i < N; ++i) {
			dual_x[i].b = 1;
			res += f(dual_x).Dual2();
			dual_x[i].b = 0;
		}

		return res;
	}

}

#endif
