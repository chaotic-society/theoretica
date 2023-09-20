
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	prec::setup("statistics");


		// Distributions

		for (int i = 1; i <= 10; ++i) {

			real alpha = i;
			real beta = 1;

			// Check mean of Gamma distribution
			prec::equals(
				"gamma_dist (1)",
				integral_gauss(
					[alpha, beta](real x) {
						return x * distribution::gamma_dist(x, alpha, beta);
					},
					tables::laguerre_roots_16, tables::laguerre_weights_16,
					16, [](real x) { return th::exp(x); }
				),
				(alpha / beta)
			);
		}

	prec::terminate();
}
