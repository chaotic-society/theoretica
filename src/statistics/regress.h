
///
/// @file regression.h Regression from data to a model.
///

#ifndef THEORETICA_REGRESSION_H
#define THEORETICA_REGRESSION_H

#include "./statistics.h"


namespace theoretica {

	/// @namespace theoretica::regression Regression to a model
	namespace regression {

		// Linear regression using least squares

		/// @class linear_model Linear regression structure
		/// for storage of least squares linear regression
		/// results with model \f$y = A + Bx\f$.
		struct linear_model {

			/// Intercept
			real A;

			/// Estimated error on A
			real err_A;

			/// Slope
			real B;

			/// Estimated error on B
			real err_B;

			/// Total error on linearization
			real err;

			/// Chi-squared on linearization
			real chi_squared;

			/// Number of degrees of freedom
			/// of the linear regression
			unsigned int ndf;

			/// The p-value associated to the
			/// computed Chi-squared
			real p_value;

		};


		/// Compute the linear regression of two sets of data of the same size
		/// using least squares linear regression. The result is a linear_regress
		/// structure with the resulting coefficients and errors. Without
		/// the error on the y axis, the chi-squared and the error
		/// on the coefficients cannot be computed.
		///
		/// @param X The set of values on the x axis
		/// @param Y The set of values on the y axis
		/// @return A linear_regress structure holding the results
		/// of the linearization
		inline linear_regress linearize(const vec_buff& X, const vec_buff& Y) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("linearize", X.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			if(Y.size() < 2) {
				TH_MATH_ERROR("linearize", Y.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			linear_regress res;

			res.A = least_squares_linear_slope(X, Y);
			res.err_A = nan();

			res.B = least_squares_linear_intercept(X, Y);
			res.err_B = nan();

			res.err = least_squares_linear_error(X, Y, res.A, res.B);
			res.chi_squared = nan();
			res.p_value = nan();
			res.ndf = Y.size() - 2;

			return res;
		}


		/// Compute the linear regression of two sets of data of the same size
		/// using least squares linear regression. The result is a linear_regress
		/// structure with the resulting coefficients and errors.
		///
		/// @param X The set of values on the x axis
		/// @param Y The set of values on the y axis
		/// @param sigma_Y The constant error on the y axis
		/// @return A linear_regress structure holding the results
		/// of the linearization
		inline linear_regress linearize(
			const vec_buff& X, const vec_buff& Y, real sigma_Y) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("linearize", X.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			if(Y.size() < 2) {
				TH_MATH_ERROR("linearize", Y.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			linear_regress res;

			res.A = least_squares_linear_slope(X, Y);
			res.err_A = least_squares_linear_sigma_A(X, Y, sigma_Y);

			res.B = least_squares_linear_intercept(X, Y);
			res.err_B = least_squares_linear_sigma_B(X, Y, sigma_Y);

			res.err = least_squares_linear_error(X, Y, res.A, res.B);
			res.chi_squared = res.err / sigma_Y;
			res.ndf = Y.size() - 2;
			res.p_value = pvalue_chi_squared(res.chi_squared, res.ndf);

			return res;
		}


		/// Compute the linear regression of two sets of data of the same size
		/// using least squares linear regression. The result is a linear_regress
		/// structure with the resulting coefficients and errors.
		///
		/// @param X The set of values on the x axis
		/// @param Y The set of values on the y axis
		/// @param sigma The different errors on the y axis
		/// @return A linear_regress structure holding the results
		/// of the linearization
		inline linear_regress linearize(
			const vec_buff& X, const vec_buff& Y, const vec_buff& sigma) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("linearize", X.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			if(Y.size() < 2) {
				TH_MATH_ERROR("linearize", Y.size(), INVALID_ARGUMENT);
				linear_regress r; r.A = nan(); r.B = nan();
				return r;
			}

			vec_buff W = sigma;
			apply(W, [](real x) {
				return 1 / square(x);
			});

			linear_regress res;

			res.A = least_squares_weighted_linear_slope(X, Y, W);
			res.err_A = 0;

			res.B = least_squares_weighted_linear_intercept(X, Y, W);
			res.err_B = 0;

			res.err = least_squares_linear_error(X, Y, res.A, res.B);
			res.chi_squared = chi_square_linearization(X, Y, sigma, res.A, res.B);
			res.ndf = Y.size() - 2;
			res.p_value = pvalue_chi_squared(res.chi_squared, res.ndf);

			return res;
		}

	}

}


#endif
