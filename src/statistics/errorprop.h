
///
/// @file erroprop.h Automatic propagation of uncertainties on arbitrary functions
///

#ifndef THEORETICA_ERRORPROP_H
#define THEORETICA_ERRORPROP_H

#include "../core/vec_buff.h"
#include "../algebra/vec.h"
#include "../autodiff/autodiff.h"
#include "../statistics/statistics.h"


namespace theoretica {


	/// Build the covariance matrix given a vector of datasets
	/// by computing the covariance between all couples of sets.
	/// @param v A vector of datasets of measures
	/// @return The covariance matrix of the datasets
	template<unsigned int N>
	inline mat<N, N> covar_mat(const std::vector<vec_buff>& v) {

		if(v.size() != N) {
			TH_MATH_ERROR("covar_mat", v.size(), INVALID_ARGUMENT);
			return mat<N, N>(nan());
		}

		mat<N, N> cm;

		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < N; ++j)
				cm.iat(i, j) = sample_covariance(v[i], v[j]);

		return cm;
	}


	/// Automatically propagate uncertainties under quadrature
	/// on an arbitrary function given the uncertainties
	/// on the variables, the mean values of the variables
	/// and the function itself, by using automatic differentiation.
	/// This function presupposes that the correlation between
	/// different variables is zero.
	///
	/// @param f The function to propagate error on
	/// @param x Best values for the variables
	/// @param delta_x Vector of uncertainties on the variables
	/// @return The propagated error on the function
	template<unsigned int N>
	inline real propagate_err(
		multidual<N>(*f)(vec<N, multidual<N>>),
		const vec<N>& x_best, const vec<N>& delta_x) {

		real err_sqr = 0;
		const multidual<N> df = f(multidual<N>::make_argument(x_best));

		for (unsigned int i = 0; i < N; ++i)
			err_sqr += square(df.Dual().get(i) * delta_x.get(i));

		return sqrt(err_sqr);
	}


	/// Automatically propagate uncertainties under quadrature
	/// on an arbitrary function given the uncertainties
	/// on the variables, the mean values of the variables
	/// and the function itself, by using automatic differentiation.
	///
	/// @param f The function to propagate error on
	/// @param x Best values for the variables
	/// @param cm Covariance matrix of the variables,
	/// where diagonal entries are the variance of the variables
	/// and off-diagonal entries are the covariance between
	/// different variables. May be constructed from datasets
	/// using the function covar_mat.
	/// @return The propagated error on the function
	template<unsigned int N>
	inline real propagate_err(
		multidual<N>(*f)(vec<N, multidual<N>>),
		const vec<N>& x_best, const mat<N, N>& cm) {

		real err_sqr = 0;
		const multidual<N> df = f(multidual<N>::make_argument(x_best));

		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < N; ++j)
				err_sqr += df.Dual().get(i) * df.Dual().get(j)
					* cm.iget(i, j);

		return sqrt(err_sqr);
	}


	/// Automatically propagate uncertainties under quadrature
	/// on an arbitrary function given the function and the
	/// set of measured data. The covar_mat function is used
	/// to estimate the covariance matrix from the data sets.
	/// For this to work, the data sets should have the same size,
	/// so as to estimate their covariance.
	///
	/// @param f The function to propagate error on
	/// @param v A vector of different datasets of the
	/// measures of the variables
	/// @return The propagated error on the function
	template<unsigned int N>
	inline real propagate_err(
		multidual<N>(*f)(vec<N, multidual<N>>),
		const std::vector<vec_buff>& v) {

		vec<N> x_mean;

		for (unsigned int i = 0; i < N; ++i)
			x_mean[i] = mean(v[i]);

		return propagate_err(f, x_mean, covar_mat<N>(v));
	}

}

#endif
