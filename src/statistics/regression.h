
///
/// @file regression.h Regression to a model.
///

#ifndef THEORETICA_REGRESSION_H
#define THEORETICA_REGRESSION_H

#include "./statistics.h"


namespace theoretica {

	/// @namespace theoretica::regression Regression to a model
	namespace regression {


		/// Compute the coefficients of the linear regression
		/// using Ordinary Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear(
			const Dataset1& X, const Dataset2& Y,
			real& intercept, real& slope, real& error) {

			// Check that the two data sets have the same size
			if(X.size() != Y.size()) {
				TH_MATH_ERROR("regression::ols_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan();
				slope = nan();
				error = nan();
				return;
			}

			// Pre-compute values
			const real sum_sqr_x = sum_squares(X);
			const real sqr_sum_x = square(sum(X));
			const real Delta = X.size() * sum_sqr_x - sqr_sum_x;
			const real prod_sum_xy = product_sum(X, Y);
			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			// Ordinary Least Squares formula for the linear model
			// without weights.
			intercept = (sum_sqr_x * sum_y - sum_x * prod_sum_xy) / Delta;
			slope = (X.size() * prod_sum_xy - sum_x * sum_y) / Delta;

			// Compute the error on the model
			real residual = 0;
			for (unsigned int i = 0; i < X.size(); ++i)
				residual += square(Y[i] - intercept - slope * X[i]);

			// Correction by degrees of freedom (N - 2)
			error = sqrt(residual / (real) (X.size() - 2));
		}


		/// Compute the coefficients of the linear regression
		/// using Ordinary Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear(
			const Dataset1& X, const Dataset2& Y, real sigma_Y,
			real& intercept, real& slope, real& error,
			real& sigma_A, real& sigma_B) {

			// Check that the two data sets have the same size
			if(X.size() != Y.size()) {
				TH_MATH_ERROR("regression::ols_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan();
				slope = nan();
				error = nan();
				return;
			}

			// Pre-compute values
			const real sum_sqr_x = sum_squares(X);
			const real sqr_sum_x = square(sum(X));
			const real Delta = X.size() * sum_sqr_x - sqr_sum_x;
			const real prod_sum_xy = product_sum(X, Y);
			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			// Ordinary Least Squares formula for the linear model
			// without weights.
			intercept = (sum_sqr_x * sum_y - sum_x * prod_sum_xy) / Delta;
			slope = (X.size() * prod_sum_xy - sum_x * sum_y) / Delta;

			// Compute the error on the model
			real residual = 0;
			for (unsigned int i = 0; i < X.size(); ++i)
				residual += square(Y[i] - intercept - slope * X[i]);

			// Correction by degrees of freedom (N - 2)
			error = sqrt(residual / (real) (X.size() - 2));

			sigma_A = sqrt(sum_sqr_x / Delta) * abs(sigma_Y);
			sigma_B = sqrt(X.size() / Delta) * abs(sigma_Y);
		}


		/// @class linear_model Linear regression structure
		/// for computation and storage of least squares linear
		/// regression results with model \f$y = A + Bx\f$.
		struct linear_model {

			/// Intercept
			real A;

			/// Estimated error on A
			real sigma_A;

			/// Slope
			real B;

			/// Estimated error on B
			real sigma_B;

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


			/// Default constructor
			linear_model() : A(nan()), B(nan()) {}


			/// Construct a linear model from data
			/// and compute the fit
			template<typename Dataset1, typename Dataset2>
			inline linear_model(const Dataset1& X, const Dataset2& Y) {
				fit(X, Y);
			}


			/// Construct a linear model from data
			/// and compute the fit
			template<typename Dataset1, typename Dataset2>
			inline linear_model(const Dataset1& X, const Dataset2& Y, real sigma_Y) {
				fit(X, Y, sigma_Y);
			}


			/// Compute the linear regression of two sets of data of the same size
			/// using ordinary least squares linear regression. Without
			/// the error on the y axis, the chi-squared and the error
			/// on the coefficients cannot be computed.
			///
			/// @param X The set of values on the x axis
			/// @param Y The set of values on the y axis
			template<typename Dataset1, typename Dataset2>
			inline void fit(const Dataset1& X, const Dataset2& Y) {

				if(X.size() != Y.size()) {
					TH_MATH_ERROR("fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() < 2) {
					TH_MATH_ERROR("fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				ols_linear(X, Y, A, B, err);

				sigma_A = nan();
				sigma_B = nan();
				chi_squared = nan();
				p_value = nan();
				ndf = Y.size() - 2;
			}


			/// Compute the linear regression of two sets of data of the same size
			/// using ordinary least squares linear regression.
			///
			/// @param X The set of values on the x axis
			/// @param Y The set of values on the y axis
			/// @param sigma_Y The constant error on the y axis
			template<typename Dataset1, typename Dataset2>
			inline void fit(
				const Dataset1& X, const Dataset2& Y, real sigma_Y) {

				if(X.size() != Y.size()) {
					TH_MATH_ERROR("fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() < 2) {
					TH_MATH_ERROR("fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(abs(sigma_Y) <= MACH_EPSILON) {
					TH_MATH_ERROR("fit", sigma_Y, DIV_BY_ZERO);
					A = nan(); B = nan();
					return;
				}

				ols_linear(X, Y, sigma_Y, A, B, err, sigma_A, sigma_B);

				chi_squared = err / sigma_Y;
				ndf = Y.size() - 2;
				p_value = pvalue_chi_squared(chi_squared, ndf);
			}


		};


		/// Compute the linear regression of two sets of data of the same size
		/// using least squares linear regression. The result is a linear_model
		/// structure with the resulting coefficients and errors.
		///
		/// @param X The set of values on the x axis
		/// @param Y The set of values on the y axis
		/// @param sigma The different errors on the y axis
		/// @return A linear_model structure holding the results
		/// of the linearization
		// inline linear_model linearize(
		// 	const vec_buff& X, const vec_buff& Y, const vec_buff& sigma) {

		// 	if(X.size() != Y.size()) {
		// 		TH_MATH_ERROR("linearize", X.size(), INVALID_ARGUMENT);
		// 		linear_model r; r.A = nan(); r.B = nan();
		// 		return r;
		// 	}

		// 	if(Y.size() < 2) {
		// 		TH_MATH_ERROR("linearize", Y.size(), INVALID_ARGUMENT);
		// 		linear_model r; r.A = nan(); r.B = nan();
		// 		return r;
		// 	}

		// 	vec_buff W = sigma;
		// 	apply(W, [](real x) {
		// 		return 1 / square(x);
		// 	});

		// 	linear_model res;

		// 	res.A = least_squares_weighted_linear_slope(X, Y, W);
		// 	res.sigma_A = 0;

		// 	res.B = least_squares_weighted_linear_intercept(X, Y, W);
		// 	res.sigma_B = 0;

		// 	res.err = least_squares_linear_error(X, Y, res.A, res.B);
		// 	res.chi_squared = chi_square_linearization(X, Y, sigma, res.A, res.B);
		// 	res.ndf = Y.size() - 2;
		// 	res.p_value = pvalue_chi_squared(res.chi_squared, res.ndf);

		// 	return res;
		// }

	}

}


#endif
