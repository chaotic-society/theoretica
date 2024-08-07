
///
/// @file autodiff.h Differential operators using automatic differentiation
///

#ifndef THEORETICA_AUTODIFF_H
#define THEORETICA_AUTODIFF_H

#include "dual.h"
#include "dual2.h"
#include "multidual.h"
#include "../algebra/vec.h"
#include "../algebra/mat.h"
#include "../core/error.h"

#include <functional>



namespace theoretica {


	// Types for multivariate automatic differentiation

	/// Real type for multivariate automatic differentiation
	/// (read "differential real").
	template<unsigned int N = 0>
	using d_real = multidual<N>;

	/// Vector type for multivariate automatic differentiation
	/// (read "differential vector").
	template<unsigned int N = 0>
	using d_vec = vec<d_real<N>, N>;


	// Univariate automatic differentiation


	/// Compute the derivative of a function at the
	/// given point using univariate automatic differentiation.
	///
	/// @param f The function to differentiate,
	/// with dual argument and return value.
	/// @param x The coordinate to compute the derivative at.
	/// @return The derivative of f at x.
	inline real deriv(dual(*f)(dual), real x) {
		return f(dual(x, 1)).Dual();
	}


	/// Compute the derivative of a function
	/// using univariate automatic differentiation.
	///
	/// @param f The function to differentiate,
	/// with dual argument and return value.
	/// @return A lambda function which computes the
	/// derivative of f using automatic differentiation.
	inline auto deriv(dual(*f)(dual)) {

		return [f](real x) {
			return deriv(f, x);
		};
	}


	/// Compute the second derivative of a function at the
	/// given point using univariate automatic differentiation.
	///
	/// @param f The function to differentiate,
	/// with dual2 argument and return value.
	/// @param x The coordinate to compute the derivative at.
	/// @return The derivative of f at x.
	inline real deriv2(dual2(*f)(dual2), real x) {
		return f(dual2(x, 1, 0)).Dual2();
	}


	/// Compute the second derivative of a function
	/// using univariate automatic differentiation.
	///
	/// @param f The function to differentiate,
	/// with dual2 argument and return value.
	/// @return A lambda function which computes the
	/// derivative of f using automatic differentiation.
	inline auto deriv2(dual2(*f)(dual2)) {

		return [f](real x) {
			return deriv2(f, x);
		};
	}


	// Differential operators

	/// Compute the gradient for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation.
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a multidual number as output.
	/// @param x The point to compute the gradient at.
	/// @return The gradient of f computed at x.
	template<unsigned int N = 0>
	inline vec<real, N> gradient(
		multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

		return f(multidual<N>::make_argument(x)).Dual();
	}


	/// Get a lambda function which computes the gradient of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a multidual number as output.
	/// @return A lambda function which computes the gradient of f.
	template<unsigned int N = 0>
	inline auto gradient(multidual<N>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {
			return gradient(f, x);
		};
	}


	/// Compute the divergence for a given \f$\vec x\f$ of a function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// using automatic differentiation
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a vector of multidual numbers as output.
	/// @param x The point to compute the divergence at.
	/// @return The divergence of f at x.
	template<unsigned int N = 0>
	inline real divergence(
		multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

		multidual<N> d = f(multidual<N>::make_argument(x));

		real div = 0;
		for (unsigned int i = 0; i < d.v.size(); ++i)
			div += d.v[i];

		return div;
	}


	/// Get a lambda function which computes the divergence of a given function
	/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given \f$\vec x\f$ using automatic differentiation.
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a vector of multidual numbers as output.
	/// @return A lambda function which computes the divergence of f at x.
	template<unsigned int N = 0>
	inline auto divergence(multidual<N>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {
			return divergence(f, x);
		};
	}


	/// Compute the jacobian of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$.
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a vector of multidual numbers as output.
	/// @param x The point to compute the Jacobian at.
	/// @return The Jacobian matrix of f at x.
	template<unsigned int N = 0, unsigned int M = 0>
	inline mat<real, M, N> jacobian(
		vec<multidual<N>, M>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

		vec<multidual<N>, M> res = f(multidual<N>::make_argument(x));
		
		// Construct the jacobian matrix
		mat<real, M, N> J;
		J.resize(res.size(), x.size());

		for (unsigned int j = 0; j < J.rows(); ++j)
			for (unsigned int i = 0; i < res[j].v.size(); ++i)
				J(j, i) = res[j].v[i];

		return J;
	}


	/// Get a lambda function which computes the jacobian of a generic
	/// function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$
	/// for a given $\vec x$.
	///
	/// @param f A function with a vector of multidual numbers as input
	/// and a vector of multidual numbers as output.
	/// @return A lambda function which computes the Jacobian matrix of f.
	template<unsigned int N = 0, unsigned int M = 0>
	inline auto jacobian(
		vec<multidual<N>, M>(*f)(vec<multidual<N>, N>)) {

		return [f](vec<real, N> x) {
			return jacobian(f, x);
		};
	}


	/// Compute the curl for a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	///
	/// @param f A function with a vector of multidual numbers as
	/// input and a vector of multidual numbers as output.
	/// @param x The point to compute the curl at.
	/// @return The curl of f at x.
	template<unsigned int N = 0>
	inline vec<real, N> curl(
		vec<multidual<N>, N>(*f)(vec<multidual<N>, N>), const vec<real, N>& x) {

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


	/// Get a lambda function which computes the curl for
	/// a given \f$\vec x\f$ of a vector field
	/// defined by \f$f: \mathbb{R}^3 \rightarrow \mathbb{R}^3\f$
	/// using automatic differentiation.
	///
	/// @param f A function with a vector of multidual numbers as
	/// input and a vector of multidual numbers as output.
	/// @return A lambda function which computes the curl of f.
	template<unsigned int N = 0>
	inline auto curl(vec<multidual<N>, N>(*f)(vec<multidual<N>, N>)) {

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
	///
	/// @param f A function with a vector of multidual numbers as
	/// input and a vector of multidual numbers as output.
	/// @param x The point to compute the directional derivative at.
	/// @param v The direction of the derivative.
	/// @return The directional derivative over v of f at x.
	template<unsigned int N = 0>
	inline vec<real, N> directional_derivative(multidual<N>(*f)(vec<multidual<N>, N>),
		const vec<real, N>& x, const vec<real, N>& v) {

		return v * dot(v, gradient(f, x));
	}


	/// Get a lambda function which computes the directional derivative
	/// of a generic function of the form
	/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// @param f The function to partially differentiate
	/// @param x The point to compute the derivative at
	/// @param v The direction to compute the derivative on
	///
	/// \note In most applications, the vector v should be a unit vector,
	/// but the function does not control whether the vector has
	/// unit length or not.
	///
	/// @param f A function with a vector of multidual numbers as
	/// input and a vector of multidual numbers as output.
	/// @param v The direction of the derivative.
	/// @return A lambda function which computes the directional
	/// derivative over v of f.
	template<unsigned int N = 0>
	inline auto directional_derivative(
		multidual<N>(*f)(vec<multidual<N>, N>), const vec<real, N>& v) {

		return [f, v](vec<real, N> x) {
			return directional_derivative(f, x, v);
		};
	}


	/// Compute the Laplacian differential operator for a generic
	/// function of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given $\vec x$.
	///
	/// @param f A function taking a vector of dual2 numbers and
	/// returning a dual2 number.
	/// @param x The point to compute the Laplacian at.
	/// @return The Laplacian of f at x.
	template<unsigned int N = 0>
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


	/// Get a lambda function which computes the Laplacian
	/// differential operator for a generic function of the
	/// form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
	/// at a given $\vec x$.
	///
	/// @param f A function taking a vector of dual2 numbers and
	/// returning a dual2 number.
	/// @return A lambda function which computes the Laplacian of f at x.
	template<unsigned int N = 0>
	inline auto laplacian(dual2(*f)(vec<dual2, N>)) {

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
	/// @param f A function taking a vector of multidual numbers and
	/// returning a multidual number.
	/// @param H The Hamiltonian of the system.
	/// @param eta A vector containing N = 2K elements, where the first
	/// K elements are the coordinates and the last K elements are the
	/// conjugate momenta.
	template<unsigned int N = 0>
	inline real sturm_liouville(
		multidual<N>(*f)(vec<multidual<N>, N>),
		multidual<N>(*H)(vec<multidual<N>, N>),
		vec<real, N> eta) {

		return gradient(f, eta)
			 * mat<real, N, N>::symplectic(eta.size(), eta.size())
			 * gradient(H, eta);
	}

}

#endif
