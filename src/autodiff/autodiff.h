
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
#include "../core/core_traits.h"
#include "./autodiff_types.h"

#include <functional>


namespace theoretica {

	/// @namespace theoretica::autodiff Differential operators with automatic differentiation.
	namespace autodiff {


		// Univariate automatic differentiation


		/// Compute the derivative of a function at the
		/// given point using univariate automatic differentiation.
		///
		/// @param f The function to differentiate,
		/// with dual argument and return value.
		/// @param x The coordinate to compute the derivative at.
		/// @return The derivative of f at x.
		template <
			typename DualFunction = std::function<dual(dual)>,
			enable_dual_func<DualFunction> = true
		>
		inline real deriv(DualFunction f, real x) {
			return f(dual(x, 1.0)).Dual();
		}


		/// Get a lambda function which computes the derivative
		/// of the given function at the given point,
		/// using automatic differentiation.
		///
		/// @param f The function to differentiate,
		/// with dual argument and return value.
		/// @return A lambda function which computes the
		/// derivative of f using automatic differentiation.
		template <
			typename DualFunction = std::function<dual(dual)>,
			enable_dual_func<DualFunction> = true
		>
		inline auto deriv(DualFunction f) {

			return [f](real x) -> real {
				return autodiff::deriv(f, x);
			};
		}


		/// Compute the second derivative of a function at the
		/// given point using univariate automatic differentiation.
		///
		/// @param f The function to differentiate,
		/// with dual2 argument and return value.
		/// @param x The coordinate to compute the derivative at.
		/// @return The derivative of f at x.
		template <
			typename Dual2Function = std::function<dual2(dual2)>,
			enable_dual2_func<Dual2Function> = true
		>
		inline real deriv2(Dual2Function f, real x) {
			return f(dual2(x, 1.0, 0.0)).Dual2();
		}


		/// Get a lambda function which computes the second
		/// derivative of the given function at the given point,
		/// using automatic differentiation.
		///
		/// @param f The function to differentiate,
		/// with dual2 argument and return value.
		/// @return A lambda function which computes the
		/// derivative of f using automatic differentiation.
		template <
			typename Dual2Function = std::function<dual2(dual2)>,
			enable_dual2_func<Dual2Function> = true
		>
		inline auto deriv2(Dual2Function f) {

			return [f](real x) -> real {
				return autodiff::deriv2(f, x);
			};
		}


		// Differential operators


		/// Compute the gradient
		/// \f$\nabla f = \sum_i^n \vec e_i \frac{\partial}{\partial x_i} f(\vec x)\f$
		/// for a given \f$\vec x\f$ of a scalar field of the form
		/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
		/// using automatic differentiation. The argument function may be
		/// a function pointer or lambda function, with the type dvec_t
		/// as first argument and dreal_t return type.
		///
		/// @param f A function with a vector of multidual numbers as input
		/// and a multidual number as output.
		/// @param x The point to compute the gradient at.
		/// @return The gradient of f computed at x.
		template <
			typename Function, typename Vector = vec<real>,
			enable_scalar_field<Function> = true,
			enable_vector<Vector> = true
		>
		inline auto gradient(Function f, const Vector& x) {

			// Extract the return type of the function
			using MultidualType = return_type_t<Function>;

			constexpr size_t N = MultidualType::size_argument;
			vec<MultidualType, N> arg;
			arg.resize(x.size());

			// Iterate over each element, setting the dual part to
			// the i-th element of the canonical base
			for (unsigned int i = 0; i < x.size(); ++i) {
				arg[i] = MultidualType(
					x[i], vec<real, N>::euclidean_base(i, x.size())
				);
			}

			return f(arg).Dual();
		}


		/// Get a lambda function which computes the gradient
		/// \f$\nabla f = \sum_i^n \vec e_i \frac{\partial}{\partial x_i} f(\vec x)\f$
		/// of a given scalar field of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
		/// at \f$\vec x\f$ using automatic differentiation.
		/// The returned lambda function accepts a vec<real, N> argument.
		/// The argument function may be a function pointer or lambda function,
		/// with the type dvec_t as first argument and dreal_t return type.
		///
		/// @param f A function with a vector of multidual numbers as input
		/// and a multidual number as output.
		/// @return A lambda function which computes the gradient of f.
		template <
			typename Function,
			enable_scalar_field<Function> = true
		>
		inline auto gradient(Function f) {

			constexpr size_t N = return_type_t<Function>::size_argument;

			return [f](const vec<real, N>& x) {
				return gradient(f, x);
			};
		}


		/// Compute the divergence
		/// \f$\sum_i^n \frac{\partial}{\partial x_i} f(\vec x)\f$
		/// for a given \f$\vec x\f$ of a scalar field of the form
		/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$ using automatic differentiation.
		/// The argument function may be a function pointer or lambda function,
		/// with the type dvec_t as first argument and dreal_t return type.
		///
		/// @param f A function with a vector of multidual numbers as input
		/// and a vector of multidual numbers as output.
		/// @param x The point to compute the divergence at.
		/// @return The divergence of f at x.
		template <
			typename Function, typename Vector = vec<real>,
			enable_scalar_field<Function> = true,
			enable_vector<Vector> = true
		>
		inline real divergence(Function f, const Vector& x) {

			using MultidualT = return_type_t<Function>;
			MultidualT d = f(MultidualT::make_argument(x));

			real div = 0;
			for (unsigned int i = 0; i < d.v.size(); ++i)
				div += d.v[i];

			return div;
		}


		/// Get a lambda function which computes the divergence of a given function
		/// of the form \f$f: \mathbb{R}^N \rightarrow \mathbb{R}\f$
		/// at a given \f$\vec x\f$ using automatic differentiation.
		/// The returned lambda function accepts a vec<real, N> argument.
		/// The argument function may be a function pointer or lambda function,
		/// with the type dvec_t as first argument and dreal_t return type.
		///
		/// @param f A function with a vector of multidual numbers as input
		/// and a vector of multidual numbers as output.
		/// @return A lambda function which computes the divergence of f at x.
		template <
			typename Function,
			enable_scalar_field<Function> = true
		>
		inline auto divergence(Function f) {

			constexpr size_t N = return_type_t<Function>::size_argument;
			using Vector = vec<real, N>;

			return [f](const Vector& x) {
				return divergence(f, x);
			};
		}


		/// Compute the jacobian of a vector field of the form
		/// \f$f: \mathbb{R}^N \rightarrow \mathbb{R}^M\f$.
		///
		/// @param f A function with a vector of multidual numbers as input
		/// and a vector of multidual numbers as output.
		/// @param x The point to compute the Jacobian at.
		/// @return The Jacobian matrix of f at x.
		template <
			typename MultidualFunction, typename Vector,
			enable_vector<Vector> = true,
			enable_vector_field<MultidualFunction> = true
		>
		inline auto jacobian(MultidualFunction f, const Vector& x) {

			constexpr size_t M = return_type_t<MultidualFunction>::size_argument;
			constexpr size_t N = _internal::func_helper<MultidualFunction>::first_arg_type::size_argument;

			const vec<multidual<N>, N> res = f(multidual<N>::make_argument(x));
			
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
		template<
			typename MultidualFunction,
			enable_vector_field<MultidualFunction> = true
		>
		inline auto jacobian(MultidualFunction f) {

			constexpr size_t N = return_type_t<MultidualFunction>::size_argument;
			using Vector = vec<real, N>;

			return [f](const Vector& x) {
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
		template <
			typename MultidualFunction, typename Vector,
			enable_vector<Vector> = true,
			enable_vector_field<MultidualFunction> = true
		>
		inline auto curl(MultidualFunction f, const Vector& x) {

			Vector res;
			res.resize(3);

			if(x.size() != 3) {
				TH_MATH_ERROR("th::curl", x.size(), MathError::InvalidArgument);
				return algebra::vec_error(res);
			}

			constexpr size_t N = return_type_t<MultidualFunction>::size_argument;
			const mat<real, N, N> J = jacobian(f, x);

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
		template<
			typename MultidualFunction,
			enable_vector_field<MultidualFunction> = true
		>
		inline auto curl(MultidualFunction f) {

			constexpr size_t N = return_type_t<MultidualFunction>::size_argument;
			using Vector = vec<real, N>;

			return [f](const Vector& x) {
				return curl(f, x);
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
		template <
			typename Dual2Function, typename Vector,
			enable_vector<Vector> = true
		>
		inline real laplacian(Dual2Function f, const Vector& x) {

			// Extract size template argument of the input vector
			constexpr size_t N = _internal::func_helper<Dual2Function>
				::first_arg_type::size_argument;

			vec<dual2, N> d;
			d.resize(x.size());

			for (unsigned int i = 0; i < d.size(); ++i)
				d[i].a = x[i];

			// Accumulate second derivatives
			real res = 0;
			for (unsigned int i = 0; i < d.size(); ++i) {
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
		template <
			typename Dual2Function,
			enable_scalar_field<Dual2Function> = true
		>
		inline auto laplacian(Dual2Function f) {

			using Vector = vec<real, return_type_t<Dual2Function>::size_argument>;

			return [f](const Vector& x) {
				return laplacian(f, x);
			};
		}

	}
}

#endif
