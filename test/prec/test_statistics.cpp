
#include "theoretica.h"
#include <cmath>
#include <ctime>
#include "prec.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	prec::setup("statistics");

		output::state.outputFiles = { "test/prec/prec_statistics.csv" };
		output::state.fieldOptions["funcName"].columnWidth = 22;

		PRNG g = PRNG::xoshiro(time(nullptr));

		// Distributions

		for (int i = 1; i <= 10; ++i) {

			real alpha = i;
			real beta = 1;

			// Check mean of Gamma distribution
			prec::equals(
				"gamma (1)",
				integral_gauss(
					[alpha, beta](real x) {
						return x * distribution::gamma(x, alpha, beta);
					},
					tables::laguerre_roots_16, tables::laguerre_weights_16,
					16, [](real x) { return th::exp(x); }
				),
				(alpha / beta)
			);
		}


		// P-value

		// Error bounds are 10^-6
		const real tol = 1E-06;

		for (int i = 0; i < 10; ++i) {

			unsigned int chi = g() % 500 + 1;
			unsigned int ndf = g() % 500 + 1;

			std::stringstream str;
			str << "pvalue(" << chi << "," << ndf << ") < 1";

			prec::equals(str.str(),
				(stats::pvalue_chi_squared(chi, ndf) - 1) < tol, 1);
		}

		for (int i = 0; i < 10; ++i) {

			unsigned int ndf = g() % 500 + 1;

			std::stringstream str;
			str << "pvalue(0," << ndf << ")";

			prec::equals(str.str(),
				stats::pvalue_chi_squared(0, ndf), 1, tol);
		}

	prec::terminate();
}
