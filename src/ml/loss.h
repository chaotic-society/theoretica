
///
/// @file loss.h Loss functions for machine learning.
///

#ifndef THEORETICA_LOSS_H
#define THEORETICA_LOSS_H

namespace theoretica {

	/// @namespace theoretica::ml Machine learning module
	namespace ml {


		/// Mean Squared Error loss function
		template<typename Vector>
		inline real mean_sqr_err(Vector y_exp, Vector y_pred) {

			if(y_exp.size() != y_pred.size()) {
				TH_MATH_ERROR("mean_sqr_err", y_exp.size(), INVALID_ARGUMENT);
				return nan();
			}

			real sqr_sum = 0;
			for (int i = 0; i < y_exp.size(); ++i)
				sqr_sum += square(y_exp[i] - y_pred[i]);

			return sqr_sum / y_exp.size();
		}


		/// Mean Absolute Error loss function
		template<typename Vector>
		inline real mean_abs_err(Vector y_exp, Vector y_pred) {

			if(y_exp.size() != y_pred.size()) {
				TH_MATH_ERROR("mean_abs_err", y_exp.size(), INVALID_ARGUMENT);
				return nan();
			}

			real abs_sum = 0;
			for (int i = 0; i < y_exp.size(); ++i)
				abs_sum += abs(y_exp[i] - y_pred[i]);

			return abs_sum / y_exp.size();
		}


		/// Huber loss function
		template<typename Vector>
		inline real huber_loss(Vector y_exp, Vector y_pred, real delta) {

			if(y_exp.size() != y_pred.size()) {
				TH_MATH_ERROR("huber_loss", y_exp.size(), INVALID_ARGUMENT);
				return nan();
			}

			real sum = 0;
			for (int i = 0; i < y_exp.size(); ++i) {

				const real diff = abs(y_exp[i] - y_pred[i]);
				sum += (diff <= delta) ?
					0.5 * square(diff) :
					delta * (diff - 0.5 * delta);
			}

			return sum / y_exp.size();
		}


		/// Log-Cosh loss function
		template<typename Vector>
		inline real logcosh_loss(Vector y_exp, Vector y_pred, real delta) {

			if(y_exp.size() != y_pred.size()) {
				TH_MATH_ERROR("huber_loss", y_exp.size(), INVALID_ARGUMENT);
				return nan();
			}

			real sum = 0;
			for (int i = 0; i < y_exp.size(); ++i)
				sum += ln(cosh(y_exp[i] - y_pred[i]));

			return sum / y_exp.size();
		}

	}

}

#endif
