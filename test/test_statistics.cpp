
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {


	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(10000);
	pdf_sampler unif = pdf_sampler::uniform(1, 20, g);

	prec::state.outputFolder = "test/";
	prec::setup("statistics");


		// Distributions

		{
			real alpha = std::floor(unif());
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
				(alpha / beta),
				0.1 // The integral computation has dominant error
			);
		}

	prec::terminate();
}
