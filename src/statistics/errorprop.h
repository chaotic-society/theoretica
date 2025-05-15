
///
/// @file errorprop.h Automatic propagation of uncertainties on arbitrary functions
///

#ifndef THEORETICA_ERRORPROP_H
#define THEORETICA_ERRORPROP_H

#include "../core/core_traits.h"
#include "../algebra/vec.h"
#include "../autodiff/autodiff.h"
#include "../statistics/statistics.h"
#include "../pseudorandom/montecarlo.h"
#include "../pseudorandom/sampling.h"


namespace theoretica {


	namespace stats {


		/// Build the covariance matrix given a vector of datasets
		/// by computing the covariance between all couples of sets.
		///
		/// @param v A vector of datasets of measures
		/// @return The covariance matrix of the datasets
		template <
			typename Matrix = mat<real>, typename Dataset = vec<real>,
			enable_vector<Dataset> = true
		>
		inline Matrix covar_mat(const std::vector<Dataset>& v) {

			Matrix cm;
			cm.resize(v.size(), v.size());

			for (unsigned int i = 0; i < cm.rows(); ++i)
				for (unsigned int j = 0; j < cm.cols(); ++j)
					cm(i, j) = stats::covariance(v[i], v[j]);

			return cm;
		}


		/// Automatically propagate uncertainties under quadrature
		/// on an arbitrary function given the uncertainties
		/// on the variables, the mean values of the variables
		/// and the function itself, by using automatic differentiation.
		/// This function assumes that the correlation between
		/// different variables is zero, if that is not the case, the covariance
		/// matrix should be used.
		///
		/// @param f The function to propagate error on
		/// @param x Best values for the variables
		/// @param delta_x Vector of uncertainties on the variables
		/// @return The propagated error on the function
		template<
			unsigned int N = 0,
			typename MultiDualFunction = autodiff::dreal_t<N>(*)(autodiff::dvec_t<N>)>
		inline real propagerr(
			MultiDualFunction f,
			const vec<real, N>& x_best, const vec<real, N>& delta_x) {

			real err_sqr = 0;
			const multidual<N> df = f(multidual<N>::make_argument(x_best));

			for (unsigned int i = 0; i < x_best.size(); ++i)
				err_sqr += square(df.Dual(i) * delta_x[i]);

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
		template <
			unsigned int N = 0,
			typename Matrix, enable_matrix<Matrix> = true,
			typename MultiDualFunction = autodiff::dreal_t<N>(*)(autodiff::dvec_t<N>)
		>
		inline real propagerr(
			MultiDualFunction f, const vec<real, N>& x_best, const Matrix& cm) {


			if(cm.rows() != x_best.size()) {
				TH_MATH_ERROR("propagerr", cm.rows(), INVALID_ARGUMENT);
				return nan();
			}

			if(cm.cols() != x_best.size()) {
				TH_MATH_ERROR("propagerr", cm.cols(), INVALID_ARGUMENT);
				return nan();
			}

			real err_sqr = 0;
			const multidual<N> df = f(multidual<N>::make_argument(x_best));

			for (unsigned int i = 0; i < cm.rows(); ++i)
				for (unsigned int j = 0; j < cm.cols(); ++j)
					err_sqr += df.Dual(i) * df.Dual(j) * cm(i, j);

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
		template<
			unsigned int N = 0,
			typename MultiDualFunction = multidual<N>(*)(autodiff::dvec_t<N>),
			typename Dataset = vec<real, N>
		>
		inline real propagerr(
			MultiDualFunction f,
			const std::vector<Dataset>& v) {

			vec<real, N> x_mean;
			x_mean.resize(v.size());

			for (unsigned int i = 0; i < v.size(); ++i)
				x_mean[i] = stats::mean(v[i]);

			return propagerr(
				f, x_mean, covar_mat<mat<real, N, N>, Dataset>(v));
		}


		/// Propagate the statistical error on a given function
		/// using the Monte Carlo method, by
		/// generating a sample following the probability
		/// distribution of the function and computing
		/// its standard deviation. N sample vectors of size M are generated
		/// by sampling the M different pdf_sampler distributions which
		/// correspond to the input variables of the function.
		/// The resulting sample is used to estimate the standard deviation
		/// over the result of the function.
		/// 
		///
		/// @param f The function to propagate error on
		/// @param rv A list of distribution samplers
		/// which sample from the probability distributions
		/// of the random variables.
		/// @param N The number of sampled values to use, defaults to
		/// 1 million.
		/// @return The standard deviation of the Monte Carlo sample
		template<typename Function>
		real propagerr_mc(
			Function f, std::vector<pdf_sampler>& rv, unsigned int N = 1E+6) {

			return stats::stdev(sample_mc(f, rv, N));
		}
	}
}


#endif
