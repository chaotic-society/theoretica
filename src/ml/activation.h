
///
/// @file activation.h Neuron activation functions
///

#ifndef THEORETICA_ACTIVATION_H
#define THEORETICA_ACTIVATION_H

#include "../core/real_analysis.h"

namespace theoretica {

	namespace ml {


		/// Rectified Linear Unit activation function
		real ReLU(real x) {
			return max(0, x);
		}


		/// Leaky Rectified Linear Unit activation function
		real leaky_ReLU(real x) {

			if(x > 0)
				return x;

			return 0.01 * x;
		}


		/// Parametric Rectified Linear Unit activation function
		real param_ReLU(real x, real a) {

			if(x > 0)
				return x;

			return a * x;
		}


		/// Swish neuron activation function
		real swish(real x) {
			return x * sigmoid(x);
		}


		/// Softplus neuron activation function
		real softplus(real x) {
			return ln(1 + exp(x));
		}

	}
}

#endif
