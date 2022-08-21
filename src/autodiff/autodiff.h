
///
/// @file autodiff.h Differential operators using automatic differentiation
///

#ifndef THEORETICA_AUTODIFF_H
#define THEORETICA_AUTODIFF_H

#include "dual.h"
#include "multidual.h"
#include "dual2.h"
#include "../algebra/vec.h"
#include "../algebra/mat.h"

#include <functional>


namespace theoretica {


	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	template<unsigned int N>
	inline vec<N> gradient(multidual<N>(*f)(vec<N, multidual<N>>), const vec<N, real>& x) {
		return f(multidual<N>::pack_function_arg(x)).Dual();
	}


	/// Get a function which computes the gradient of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<vec<N, real>(vec<N, real>)> gradient(multidual<N>(*f)(vec<N, multidual<N>>)) {

		return [f](vec<N, real> x) {
			return f(multidual<N>::pack_function_arg(x)).Dual();
		};
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

		for (unsigned int i = 0; i < N; ++i)
			dual_x[i] = x.get(i);

		for (unsigned int i = 0; i < N; ++i) {
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
		for (unsigned int i = 0; i < N; ++i)
			div += d.v.at(i);

		return div;
	}


	/// Get a function which computes the divergence of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<real(vec<N, real>)>
	divergence(multidual<N>(*f)(vec<N, multidual<N>>)) {

		return [f](vec<N, real> x) {

			multidual<N> d = f(multidual<N>::pack_function_arg(x));

			real div = 0;
			for (unsigned int i = 0; i < N; ++i)
				div += d.v.at(i);

			return div;
		};
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

		for (unsigned int i = 0; i < N; ++i)
			dual_x[i] = x.get(i);

		for (unsigned int i = 0; i < N; ++i) {
			dual_x[i].b = 1;
			res += f(dual_x).Dual();
			dual_x[i].b = 0;
		}

		return res;
	}


	/// Compute the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$
	template<unsigned int N, unsigned int M>
	inline mat<M, N> jacobian(vec<M, multidual<N>>(*f)(vec<N, multidual<N>>), const vec<N, real>& x) {

		vec<M, multidual<N>> res = f(multidual<N>::pack_function_arg(x));

		// Construct the jacobian matrix
		mat<M, N> J;
		for (unsigned int i = 0; i < N; ++i) {
			for (unsigned int j = 0; j < M; ++j) {
				J.iat(j, i) = res.at(j).Dual().at(i);
			}
		}

		return J;
	}


	/// Get a function which computes the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$ for a given $\vec x$
	template<unsigned int N, unsigned int M>
	inline std::function<mat<M, N>(vec<N, real>)> jacobian(vec<M, multidual<N>>(*f)(vec<N, multidual<N>>)) {

		return [f](vec<N, real> x) {

			vec<M, multidual<N>> res = f(multidual<N>::pack_function_arg(x));

			// Construct the jacobian matrix
			mat<M, N> J;
			for (unsigned int i = 0; i < N; ++i) {
				for (unsigned int j = 0; j < M; ++j) {
					J.iat(j, i) = res.at(j).Dual().at(i);
				}
			}

			return J;
		};
	}


	/// Compute the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	inline vec<3> curl(vec<3, multidual<3>>(*f)(vec<3, multidual<3>>), const vec3& x) {

		mat3 J = jacobian<3, 3>(f, x);
		vec3 res;

		res.at(0) = J.iat(2, 1) - J.iat(1, 2);
		res.at(1) = J.iat(0, 2) - J.iat(2, 0);
		res.at(2) = J.iat(1, 0) - J.iat(0, 1);

		return res;
	}


	/// Get a function which computes the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	inline std::function<vec3(vec3)> curl(vec<3, multidual<3>>(*f)(vec<3, multidual<3>>)) {

		return [f](vec3 x) {

			mat3 J = jacobian<3, 3>(f, x);
			vec3 res;

			res.at(0) = J.iat(2, 1) - J.iat(1, 2);
			res.at(1) = J.iat(0, 2) - J.iat(2, 0);
			res.at(2) = J.iat(1, 0) - J.iat(0, 1);

			return res;

		};
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

		// Remove .normalized() to avoid sqrt() ?

		return v.normalized() * dot(v, gradient(f, x));
	}


	/// Get a function which computes the directional derivative of a generic function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// @param f The function to partially differentiate
	/// @param x The point to compute the derivative at
	/// @param v The direction to compute the derivative on
	///
	/// \note In most applications, the vector v should be a unit vector,
	/// but the function does not control whether the vector has
	/// unit length or not.
	template<unsigned int N>
	inline std::function<vec<N, real>(vec<N, real>)>
	directional_derivative(multidual<N>(*f)(vec<N, multidual<N>>), const vec<N, real>& v) {

		// Remove .normalized() to avoid sqrt() ?

		return [f, v](vec<N, real> x) {
			return v.normalized() * dot(v, gradient(f, x));
		};
	}


	/// Compute the laplacian differential operator for a generic
	/// function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given $\vec x$.
	template<unsigned int N>
	inline real laplacian(dual2(*f)(vec<N, dual2>), const vec<N>& x) {

		real res = 0;
		vec<N, dual2> d;

		for (unsigned int i = 0; i < N; ++i)
			d[i] = x.get(i);

		for (unsigned int i = 0; i < N; ++i) {
			d[i].b = 1;
			res += f(d).Dual2();
			d[i].b = 0;
		}

		return res;
	}


	/// Get a function which computes the laplacian differential operator
	/// for a generic function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given $\vec x$.
	template<unsigned int N>
	inline std::function<real(vec<N, real>)> laplacian(dual2(*f)(vec<N, dual2>)) {

		return [f](vec<N, real> x) {

			real res = 0;
			vec<N, dual2> d;

			for (unsigned int i = 0; i < N; ++i)
				d[i] = x.get(i);

			for (unsigned int i = 0; i < N; ++i) {
				d[i].b = 1;
				res += f(d).Dual2();
				d[i].b = 0;
			}

			return res;
		};
	}


	/// Compute the Sturm-Liouville operator on a generic function
	/// of the form \f$f: \mathbb{R}^{2N} \rightarrow \mathbb{R}\f$
	/// with respect to a given Hamiltonian function of the form
	/// \f$H: \mathbb{R}^{2N} \rightarrow \mathbb{R}\f$ where
	/// the first N arguments are the coordinates in phase space
	/// and the last N arguments are the conjugate momenta,
	/// for a given point in phase space.
	///
	/// @param f The function to apply the operator to
	/// @param H The Hamiltonian
	/// @param eta A vector containing M = 2N elements, where the first
	/// N elements are the coordinates and the last N elements are the
	/// conjugate momenta.
	template<unsigned int M>
	inline real sturm_liouville(dual(*f)(vec<M, dual>), dual(*H)(vec<M, dual>), vec<M> eta) {
		return gradient(f, eta) * mat<M, M>::symplectic() * gradient(H, eta);
	}

}

#endif
