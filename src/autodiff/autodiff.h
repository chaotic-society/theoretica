
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

#include "../core/error.h"


namespace theoretica {


	/// Types for multivariate automatic differentiation

	template<unsigned int N = 0>
	using d_real = multidual<N>;

	template<unsigned int N = 0>
	using d_vec = vec<d_real<N>, N>;


	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	template<unsigned int N>
	inline vec<real, N> gradient(
		d_real<N>(*f)(d_vec<N>), const vec<real, N>& x) {

		return f(d_real<N>::make_argument(x)).Dual();
	}


	/// Get a function which computes the gradient of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<vec<real, N>(vec<real, N>)> gradient(
		d_real<N>(*f)(d_vec<N>)) {

		return [f](vec<real, N> x) {
			return gradient(f, x);
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
		res.resize(x.size());

		vec<dual, N> dual_x;
		dual_x.resize(x.size());

		for (unsigned int i = 0; i < x.size(); ++i)
			dual_x[i].a = x[i];

		// Re-evaluate the function with varying
		// dual part to extract each derivative
		for (unsigned int i = 0; i < x.size(); ++i) {
			dual_x[i].b = 1;
			res[i] = f(dual_x).Dual();
			dual_x[i].b = 0;
		}

		return res;
	}


	/// Compute the divergence for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation
	template<unsigned int N>
	inline real divergence(d_real<N>(*f)(d_vec<N>), const vec<real, N>& x) {

		d_real<N> d = f(d_real<N>::make_argument(x));

		real div = 0;
		for (unsigned int i = 0; i < d.v.size(); ++i)
			div += d.v[i];

		return div;
	}


	/// Get a function which computes the divergence of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	template<unsigned int N>
	inline std::function<real(vec<real, N>)>
	divergence(d_real<N>(*f)(d_vec<N>)) {

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
		dual_x.resize(x.size());

		for (unsigned int i = 0; i < N; ++i)
			dual_x[i].a = x[i];

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
		vec<d_real<N>, M>(*f)(d_vec<N>), const vec<real, N>& x) {

		vec<d_real<N>, M> res = f(d_real<N>::make_argument(x));
		
		// Construct the jacobian matrix
		mat<real, M, N> J;
		J.resize(res.size(), x.size());

		for (unsigned int j = 0; j < J.rows(); ++j)
			for (unsigned int i = 0; i < res[j].v.size(); ++i)
				J(j, i) = res[j].v[i];

		return J;
	}


	/// Get a function which computes the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$ for a given $\vec x$
	template<unsigned int N, unsigned int M>
	inline std::function<mat<real, M, N>(vec<real, N>)> jacobian(
		vec<d_real<N>, M>(*f)(d_vec<N>)) {

		return [f](vec<real, N> x) {
			return jacobian(f, x);
		};
	}


	/// Compute the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	template<unsigned int N>
	inline vec<real, N> curl(d_vec<N>(*f)(d_vec<N>), const vec<real, N>& x) {

		if(x.size() != 3) {
			TH_MATH_ERROR("th::curl", x.size(), INVALID_ARGUMENT);
			return vec<real, N>(nan(), x.size());
		}

		mat<real, N, N> J = jacobian<N, N>(f, x);
		J.resize(3, 3);

		vec<real, N> res;
		res.resize(3);

		res[0] = J(2, 1) - J(1, 2);
		res[1] = J(0, 2) - J(2, 0);
		res[2] = J(1, 0) - J(0, 1);

		return res;
	}


	/// Get a function which computes the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	template<unsigned int N>
	inline std::function<vec<real, N>(vec<real, N>)> curl(d_vec<N>(*f)(d_vec<N>)) {

		return [f](vec<real, N> x) {
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
	inline vec<real, N> directional_derivative(d_real<N>(*f)(d_vec<N>),
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
	directional_derivative(d_real<N>(*f)(d_vec<N>), const vec<real, N>& v) {

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
		d.resize(x.size());

		for (unsigned int i = 0; i < x.size(); ++i)
			d[i].a = x[i];

		for (unsigned int i = 0; i < x.size(); ++i) {
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
	/// @param eta A vector containing N = 2K elements, where the first
	/// K elements are the coordinates and the last K elements are the
	/// conjugate momenta.
	template<unsigned int N>
	inline real sturm_liouville(
		d_real<N>(*f)(d_vec<N>), d_real<N>(*H)(d_vec<N>), vec<real, N> eta) {

		return gradient(f, eta)
			 * mat<real, N, N>::symplectic(eta.size(), eta.size())
			 * gradient(H, eta);
	}

}

#endif
