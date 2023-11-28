
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
	inline vec<real, N> gradient(multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {
		return f(multidual<N>::make_argument(x)).Dual();
	}


	/// Get a function which computes the gradient of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<vec<real, N>(vec<real, N>)> gradient(multidual<N>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {
			return f(multidual<N>::make_argument(x)).Dual();
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
	inline vec<real, N> gradient_mono(dual(*f)(vec<dual, N>), const vec<real, N>& x) {

		vec<real, N> res;
		vec<dual, N> dual_x;

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
	inline real divergence(multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

		multidual<N> d = f(multidual<N>::make_argument(x));

		real div = 0;
		for (unsigned int i = 0; i < N; ++i)
			div += d.v.get(i);

		return div;
	}


	/// Get a function which computes the divergence of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<real(vec<real, N>)>
	divergence(multidual<N>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {

			return divergence(f, x);
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
	inline real divergence_mono(dual(*f)(vec<dual, N>), const vec<real, N>& x) {

		real res = 0;
		vec<dual, N> dual_x;

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
	inline mat<real, M, N> jacobian(
		vec<multidual<N>, M>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

		vec<multidual<N>, M> res = f(multidual<N>::make_argument(x));

		// Construct the jacobian matrix
		mat<real, M, N> J;
		for (unsigned int i = 0; i < N; ++i) {
			for (unsigned int j = 0; j < M; ++j) {
				J.at(j, i) = res.at(j).Dual().at(i);
			}
		}

		return J;
	}


	/// Get a function which computes the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$ for a given $\vec x$
	template<unsigned int N, unsigned int M>
	inline std::function<mat<real, M, N>(vec<real, N>)> jacobian(
		vec<multidual<N>, M>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {
			return jacobian(f, x);
		};
	}


	/// Compute the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	inline vec3 curl(vec<multidual<3>, 3>(*f)(vec<multidual<3>, 3>), const vec3& x) {

		mat3 J = jacobian<3, 3>(f, x);
		vec3 res;

		res.at(0) = J.at(2, 1) - J.at(1, 2);
		res.at(1) = J.at(0, 2) - J.at(2, 0);
		res.at(2) = J.at(1, 0) - J.at(0, 1);

		return res;
	}


	/// Get a function which computes the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	inline std::function<vec3(vec3)> curl(vec<multidual<3>, 3>(*f)(vec<multidual<3>, 3>)) {

		return [f](vec3 x) {
			return curl(f, x);
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
	inline vec<real, N> directional_derivative(multidual<N>(*f)(vec<multidual<N>, N>),
		const vec<real, N>& x, const vec<real, N>& v) {

		return v * dot(v, gradient(f, x));
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
	inline std::function<vec<real, N>(vec<real, N>)>
	directional_derivative(multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& v) {

		return [f, v](vec<real, N> x) {
			return directional_derivative(f, x, v);
		};
	}


	/// Compute the laplacian differential operator for a generic
	/// function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given $\vec x$.
	template<unsigned int N>
	inline real laplacian(dual2(*f)(vec<dual2, N>), const vec<real, N>& x) {

		real res = 0;
		vec<dual2, N> d;

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
	inline std::function<real(vec<real, N>)> laplacian(dual2(*f)(vec<dual2, N>)) {

		return [f](vec<real, N> x) {
			return laplacian(f, x);
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
	inline real sturm_liouville(dual(*f)(vec<dual, M>), dual(*H)(vec<dual, M>), vec<real, M> eta) {
		return gradient(f, eta) * mat<real, M, M>::symplectic() * gradient(H, eta);
	}

}

#endif
