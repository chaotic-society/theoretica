
///
/// @file erroprop.h Automatic propagation of uncertainties on arbitrary functions
///

#ifndef THEORETICA_ERRORPROP_H
#define THEORETICA_ERRORPROP_H

#include "../core/vec_buff.h"
#include "../algebra/vec.h"
#include "../autodiff/autodiff.h"
#include "../statistics/statistics.h"
#include "../pseudorandom/montecarlo.h"
#include "../pseudorandom/rand_dist.h"


namespace theoretica {


	/// Build the covariance matrix given a vector of datasets
	/// by computing the covariance between all couples of sets.
	/// @param v A vector of datasets of measures
	/// @return The covariance matrix of the datasets
	template<unsigned int N>
	inline mat<real, N, N> covar_mat(const std::vector<vec_buff>& v) {

		if(v.size() != N) {
			TH_MATH_ERROR("covar_mat", v.size(), INVALID_ARGUMENT);
			return mat<real, N, N>(nan());
		}

		mat<real, N, N> cm;

		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < N; ++j)
				cm.at(i, j) = sample_covariance(v[i], v[j]);

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
	inline real error_propagation(
		multidual<N>(*f)(vec<multidual<N>, N>),
		const vec<real, N>& x_best, const vec<real, N>& delta_x) {

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
	inline real error_propagation(
		multidual<N>(*f)(vec<multidual<N>, N>),
		const vec<real, N>& x_best, const mat<real, N, N>& cm) {

		real err_sqr = 0;
		const multidual<N> df = f(multidual<N>::make_argument(x_best));

		for (unsigned int i = 0; i < N; ++i)
			for (unsigned int j = 0; j < N; ++j)
				err_sqr += df.Dual().get(i) * df.Dual().get(j)
					* cm.get(i, j);

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
	inline real error_propagation(
		multidual<N>(*f)(vec<multidual<N>, N>),
		const std::vector<vec_buff>& v) {

		vec<real, N> x_mean;

		for (unsigned int i = 0; i < N; ++i)
			x_mean[i] = mean(v[i]);

		return error_propagation(f, x_mean, covar_mat<N>(v));
	}


	/// Propagate the statistical error on a given function
	/// using the Monte Carlo method, by
	/// generating a sample following the probability
	/// distribution of the function and computing
	/// its standard deviation.
	///
	/// @param f The function to propagate error on
	/// @param rv A list of distribution samplers
	/// which sample from the probability distributions
	/// of the random variables.
	/// @param N The number of sampled values to use, defaults to
	/// 1 million.
	/// @return The standard deviation of the Monte Carlo sample
	template<typename Function>
	real mc_error_propagation(
		Function f, std::vector<pdf_sampler> rv, unsigned int N = 1E+6) {

		return smpl_stdev(mc_sample(f, rv, N));
	}

}

#endif
