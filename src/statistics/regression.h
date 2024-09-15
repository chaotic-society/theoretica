
///
/// @file regression.h Regression to a model.
///

#ifndef THEORETICA_REGRESSION_H
#define THEORETICA_REGRESSION_H

#include "./statistics.h"
#include "../algebra/mat.h"


namespace theoretica {

	/// @namespace theoretica::regression Regression to a model
	namespace regression {


		/// Compute the coefficients of the linear regression
		/// using Ordinary Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear(
			const Dataset1& X, const Dataset2& Y,
			real& intercept, real& slope) {

			// Check that the two data sets have the same size
			if(X.size() != Y.size()) {
				TH_MATH_ERROR("regression::ols_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan(); slope = nan();
				return;
			}

			// Pre-compute values
			const real sum_sqr_x = sum_squares(X);
			const real sqr_sum_x = square(sum(X));
			const real Delta = X.size() * sum_sqr_x - sqr_sum_x;
			const real prod_sum_xy = product_sum(X, Y);
			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			if(abs(Delta) < MACH_EPSILON) {
				TH_MATH_ERROR("ols_linear", Delta, DIV_BY_ZERO);
				intercept = nan(); slope = nan();
				return;
			}

			// Ordinary Least Squares formula for the linear model
			// without weights.
			intercept = (sum_sqr_x * sum_y - sum_x * prod_sum_xy) / Delta;
			slope = (X.size() * prod_sum_xy - sum_x * sum_y) / Delta;
		}


		/// Compute the coefficients of the linear regression
		/// using Ordinary Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear(
			const Dataset1& X, const Dataset2& Y, real sigma_Y,
			real& intercept, real& slope,
			real& sigma_A, real& sigma_B) {

			// Check that the two data sets have the same size
			if(X.size() != Y.size()) {
				TH_MATH_ERROR("regression::ols_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan(); slope = nan();
				return;
			}

			// Pre-compute values
			const real sum_sqr_x = sum_squares(X);
			const real sqr_sum_x = square(sum(X));
			const real Delta = X.size() * sum_sqr_x - sqr_sum_x;
			const real sum_xy = product_sum(X, Y);
			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			if(abs(Delta) < MACH_EPSILON) {
				TH_MATH_ERROR("ols_linear", Delta, DIV_BY_ZERO);
				intercept = nan(); slope = nan();
				return;
			}

			// Ordinary Least Squares formula for the linear model
			// without weights.
			intercept = (sum_sqr_x * sum_y - sum_x * sum_xy) / Delta;
			slope = (X.size() * sum_xy - sum_x * sum_y) / Delta;

			// Compute the uncertainties on the coefficients
			sigma_A = sqrt(square(sigma_Y) * sum_sqr_x / Delta);
			sigma_B = sqrt(square(sigma_Y) * X.size() / Delta);
		}


		/// Compute the coefficients of the linear regression
		/// using Ordinary Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear(
			const Dataset1& X, const Dataset2& Y, real sigma_Y,
			real& intercept, real& slope, mat2& covar_mat) {

			// Check that the two data sets have the same size
			if(X.size() != Y.size()) {
				TH_MATH_ERROR("regression::ols_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan();
				slope = nan();
				return;
			}

			// Pre-compute values
			const real sum_sqr_x = sum_squares(X);
			const real sqr_sum_x = square(sum(X));
			const real Delta = X.size() * sum_sqr_x - sqr_sum_x;
			const real sum_xy = product_sum(X, Y);
			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			if(abs(Delta) < MACH_EPSILON) {
				TH_MATH_ERROR("ols_linear", Delta, DIV_BY_ZERO);
				intercept = nan(); slope = nan();
				return;
			}

			// Ordinary Least Squares formula for the linear model
			// without weights.
			intercept = (sum_sqr_x * sum_y - sum_x * sum_xy) / Delta;
			slope = (X.size() * sum_xy - sum_x * sum_y) / Delta;

			// Compute the Covariance Matrix for the coefficients
			covar_mat(0, 0) = square(sigma_Y) * sum_sqr_x / Delta;
			covar_mat(1, 1) = square(sigma_Y) * X.size() / Delta;
			covar_mat(0, 1) = -square(sigma_Y) * sum_x / Delta;
			covar_mat(1, 0) = covar_mat(0, 1);
		}


		/// Compute the coefficients of the linear regression
		/// using Weighted Least Squares
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline void wls_linear(
			const Dataset1& X, const Dataset2& Y, const Dataset3& W,
			real& intercept, real& slope, mat2& covar_mat) {

			if(X.size() != Y.size() || X.size() != W.size()) {
				TH_MATH_ERROR("wls_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan(); slope = nan();
				return;
			}

			// Pre-compute values
			const real sum_xw = product_sum(X, W);
			const real sum_xxw = product_sum(X, X, W);
			const real sum_xyw = product_sum(X, Y, W);
			const real sum_yw = product_sum(Y, W);
			const real sum_w = sum(W);
			const real Delta = sum_w * sum_xxw - square(sum_xw);

			if(abs(Delta) < MACH_EPSILON) {
				TH_MATH_ERROR("wls_linear", Delta, DIV_BY_ZERO);
				intercept = nan(); slope = nan();
				return;
			}

			intercept = (sum_xxw * sum_yw - sum_xw * sum_xyw) / Delta;
			slope = (sum_w * sum_xyw - sum_xw * sum_yw) / Delta;

			// Compute the Covariance Matrix for the coefficients
			covar_mat(0, 0) = sum_xxw / Delta;
			covar_mat(1, 1) = sum_w / Delta;
			covar_mat(0, 1) = -sum_xw / Delta;
			covar_mat(1, 0) = covar_mat(0, 1);
		}


		/// Compute the coefficients of the linear regression
		/// using Weighted Least Squares
		template<typename Dataset1, typename Dataset2>
		inline void wls_linear(
			const Dataset1& X, const Dataset2& Y,
			real sigma_X, real sigma_Y,
			real& intercept, real& slope,
			mat2& covar_mat) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("wls_linear", X.size(), INVALID_ARGUMENT);
				intercept = nan();
				slope = nan();
				return;
			}

			// Pre-compute values
			const real mean_x = stats::mean(X);
			const real mean_y = stats::mean(Y);
			const real delta_x = (sum_squares(X) / X.size()) - square(mean_x);
			const real delta_y = (sum_squares(Y) / Y.size()) - square(mean_y);
			const real delta_xy = (product_sum(X, Y) / X.size()) - mean_x * mean_y;
			const real Delta = square(sigma_X) * delta_y - square(sigma_Y) * delta_x;

			slope = (Delta + sqrt(
					square(Delta) + 4 * square(sigma_X * sigma_Y * delta_xy)
				)) / (2 * square(sigma_X) * delta_xy);

			intercept = mean_y - slope * mean_x;

			// TO-DO Compute Covariance Matrix for this case
			covar_mat = mat2(nan());
		}

		/// Compute the Ordinary Least Squares regression
		/// to a line passing through the origin.
		template<typename Dataset1, typename Dataset2>
		inline void ols_linear_orig(
			const Dataset1& X, const Dataset2& Y, real sigma_Y,
			real& B, real& sigma_B) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("ols_linear_orig", X.size(), INVALID_ARGUMENT);
				B = nan();
				return;
			}

			const real sum_x = sum(X);
			const real sum_y = sum(Y);

			if(abs(sum_y) < MACH_EPSILON) {
				TH_MATH_ERROR("ols_linear_orig", sum_y, DIV_BY_ZERO);
				B = nan();
				return;
			}

			B = sum_x / sum_y;
			sigma_B = sqrt(X.size() * square(sigma_Y / sum_x));
		}

		/// Compute the Weight Least Squares regression
		/// to a line passing through the origin.
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline void wls_linear_orig(
			const Dataset1& X, const Dataset2& Y,
			const Dataset3& W, real& B, real& sigma_B) {

			if(X.size() != Y.size() || X.size() != W.size()) {
				TH_MATH_ERROR("wls_linear_orig", X.size(), INVALID_ARGUMENT);
				B = nan();
				return;
			}

			const real sum_xw = product_sum(X, W);
			const real sum_yw = product_sum(Y, W);

			if(abs(sum_yw) < MACH_EPSILON) {
				TH_MATH_ERROR("ols_linear_orig", sum_yw, DIV_BY_ZERO);
				B = nan();
				return;
			}

			B = sum_xw / sum_yw;
			sigma_B = sqrt(sum(W) / sum_xw);
		}


		/// Compute the error of the least squares
		/// linear regression from the X and Y datasets
		template<typename Dataset1, typename Dataset2>
		inline real ols_linear_error(
			const Dataset1& X, const Dataset2& Y,
			real intercept, real slope) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("ols_linear_error", X.size(), INVALID_ARGUMENT);
				return nan();
			}

			real err = 0;
			for (unsigned int i = 0; i < X.size(); ++i) {
				err += square(Y[i] - intercept - slope * X[i]);
			}

			// Correction by degrees of freedom (N - 2)
			return sqrt(err / (real) (X.size() - 2));
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

			/// Covariance matrix of the coefficients A and B
			mat2 covar_mat;

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


			/// Construct a linear model from data
			/// and compute the fit
			template<typename Dataset1, typename Dataset2, typename Dataset3>
			inline linear_model(
				const Dataset1& X, const Dataset2& Y, const Dataset3& sigma_Y) {
				fit(X, Y, sigma_Y);
			}


			/// Construct a linear model from data
			/// and compute the fit
			template<typename Dataset1, typename Dataset2>
			inline linear_model(
				const Dataset1& X, const Dataset2& Y, real sigma_X, real sigma_Y) {
				fit(X, Y, sigma_X, sigma_Y);
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
					TH_MATH_ERROR("linear_model::fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() <= 2) {
					TH_MATH_ERROR("linear_model::fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				ols_linear(X, Y, A, B);

				if(is_nan(A))
					return;

				err = ols_linear_error(X, Y, A, B);
				covar_mat = mat2(nan());
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
					TH_MATH_ERROR("linear_model::fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() <= 2) {
					TH_MATH_ERROR("linear_model::fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(abs(sigma_Y) <= MACH_EPSILON) {
					TH_MATH_ERROR("linear_model::fit", sigma_Y, DIV_BY_ZERO);
					A = nan(); B = nan();
					return;
				}

				ols_linear(X, Y, sigma_Y, A, B, this->covar_mat);

				if(is_nan(A))
					return;

				sigma_A = sqrt(this->covar_mat(0, 0));
				sigma_B = sqrt(this->covar_mat(1, 1));

				err = ols_linear_error(X, Y, A, B);
				chi_squared = err / sigma_Y;
				ndf = Y.size() - 2;
				p_value = stats::pvalue_chi_squared(chi_squared, ndf);
			}


			/// Compute the linear regression of two sets of data of the same size
			/// using least squares linear regression.
			///
			/// @param X The set of values on the x axis
			/// @param Y The set of values on the y axis
			/// @param sigma The different errors on the y axis
			template<typename Dataset1, typename Dataset2, typename Dataset3>
			inline void fit(
				const Dataset1& X, const Dataset2& Y, const Dataset3& sigma) {

				if(X.size() != Y.size()) {
					TH_MATH_ERROR("linear_model::fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() < 2) {
					TH_MATH_ERROR("linear_model::fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				// The weights are given by the squared inverse of the sigma
				auto W = sigma;
				apply(W, [](real x) {
					return 1.0 / (x * x);
				});

				wls_linear(X, Y, W, A, B, covar_mat);

				if(is_nan(A))
					return;

				sigma_A = sqrt(this->covar_mat(0, 0));
				sigma_B = sqrt(this->covar_mat(1, 1));

				err = ols_linear_error(X, Y, A, B);
				chi_squared = stats::chi_square_linearization(X, Y, sigma, A, B);
				ndf = Y.size() - 2;
				p_value = stats::pvalue_chi_squared(chi_squared, ndf);
			}


			/// Compute the linear regression of two sets of data of the same size
			/// using least squares linear regression.
			///
			/// @param X The set of values on the x axis
			/// @param Y The set of values on the y axis
			/// @param sigma The different errors on the y axis
			template<typename Dataset1, typename Dataset2>
			inline void fit(
				const Dataset1& X, const Dataset2& Y, real sigma_X, real sigma_Y) {

				if(X.size() != Y.size()) {
					TH_MATH_ERROR("linear_model::fit", X.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				if(Y.size() < 2) {
					TH_MATH_ERROR("linear_model::fit", Y.size(), INVALID_ARGUMENT);
					A = nan(); B = nan();
					return;
				}

				wls_linear(X, Y, sigma_X, sigma_Y, A, B, covar_mat);

				if(is_nan(A))
					return;

				sigma_A = sqrt(this->covar_mat(0, 0));
				sigma_B = sqrt(this->covar_mat(1, 1));

				err = ols_linear_error(X, Y, A, B);
				chi_squared = err / (square(sigma_Y) + square(B * sigma_X));
				ndf = Y.size() - 2;
				p_value = stats::pvalue_chi_squared(chi_squared, ndf);
			}


			/// Compute the expected Y value for the given X
			/// value, following the computed model.
			inline real operator()(real x) {
				return A + B * x;
			}


#ifndef THEORETICA_NO_PRINT

		/// Convert the vector to string representation
		inline std::string to_string() const {

			std::stringstream res;

			res << "Model: y = A + B * x" << std::endl;

			res << "A = " << A;
			if(!is_nan(sigma_A))
				res << " +- " << sigma_A;
			res << std::endl;

			res << "B = " << B;
			if(!is_nan(sigma_B))
				res << " +- " << sigma_B;
			res << std::endl;

			res << "Error = " << err << std::endl;
			res << "ndf = " << ndf << std::endl;

			if(!is_nan(chi_squared))
				res << "Chi-squared = " << chi_squared << std::endl;

			if(!is_nan(p_value))
				res << "p-value = " << p_value << std::endl;

			if(!is_nan(covar_mat(0, 0)))
				res << "Covariance Matrix:\n" << covar_mat;

			return res.str();
		}


		/// Convert the vector to string representation.
		inline operator std::string() {
			return to_string();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(
			std::ostream& out, const linear_model& obj) {
			return out << obj.to_string();
		}

#endif

		};


		/// Compute the intercept of the least squares
		/// linear regression from X and Y
		template<typename Dataset1, typename Dataset2>
		inline real ols_linear_intercept(
			const Dataset1& X, const Dataset2& Y) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("ols_linear_intercept", X.size(), INVALID_ARGUMENT);
				return nan();
			}

			const real Delta = X.size() * sum_squares(X) - square(sum(X));
			const real A = (sum_squares(X) * sum(Y) - sum(X)
				* product_sum(X, Y)) / Delta;

			return A;
		}


		/// Compute the error on the intercept (A)
		template<typename Dataset1, typename Dataset2>
		inline real ols_linear_sigma_A(
			const Dataset1& X, const Dataset2& Y, real sigma_y) {

			const real Delta = X.size() * sum_squares(X) - square(sum(X));
			return sqrt(sum_squares(X) / Delta) * abs(sigma_y);
		}


		/// Compute the slope of the least squares
		/// linear regression from X and Y
		template<typename Dataset1, typename Dataset2>
		inline real ols_linear_slope(const Dataset1& X, const Dataset2& Y) {

			if(X.size() != Y.size()) {
				TH_MATH_ERROR("ols_linear_slope", X.size(), INVALID_ARGUMENT);
				return nan();
			}

			const real Delta = X.size() * sum_squares(X) - square(sum(X));
			const real B = (X.size() * product_sum(X, Y) - sum(X) * sum(Y))
				/ Delta;

			return B;
		}


		/// Compute the error on the slope coefficient (B)
		template<typename Dataset1, typename Dataset2>
		inline real ols_linear_sigma_B(
			const Dataset1& X, const Dataset2& Y, real sigma_y) {

			const real Delta = X.size() * sum_squares(X) - square(sum(X));
			return sqrt(X.size() / Delta) * abs(sigma_y);
		}


		/// Compute the intercept of the weighted least squares
		/// linear regression from X and Y
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline real wls_linear_intercept(
			const Dataset1& X, const Dataset2& Y, const Dataset3& W) {

			if(X.size() != Y.size() || X.size() != W.size()) {
				TH_MATH_ERROR(
					"wls_linear_intercept",
					X.size(), INVALID_ARGUMENT);
				return nan();
			}

			const real Delta = sum(W) * product_sum(X, X, W)
				- square(product_sum(X, W));

			const real A = (product_sum(X, X, W) * product_sum(Y, W) -
				product_sum(X, W) * product_sum(X, Y, W)) / Delta;

			return A;
		}


		/// Compute the slope of the weighted least squares
		/// linear regression from X and Y
		template<typename Dataset1, typename Dataset2, typename Dataset3>
		inline real wls_linear_slope(
			const Dataset1& X, const Dataset2& Y, const Dataset3& W) {

			if(X.size() != Y.size() || X.size() != W.size()) {
				TH_MATH_ERROR(
					"wls_linear_slope",
					X.size(), INVALID_ARGUMENT);
				return nan();
			}

			const real Delta = sum(W) * product_sum(X, X, W)
				- square(product_sum(X, W));

			const real B = (sum(W) * product_sum(X, Y, W) -
				product_sum(X, W) * product_sum(Y, W)) / Delta;

			return B;
		}
	}
}


#endif
