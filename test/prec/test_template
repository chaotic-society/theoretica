
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	// Start testing the module
	prec::setup("moduleName");

		// Set the output file for the tests
		output::options.outputFiles = { "test/prec/prec_module.csv" };

		/// Check that x_approx is almost equal to x_exact
		prec::equals("funcName", x_approx, x_exact, options);

		// Check that f_approx is almost equal to f_exact over a domain,
		// estimating different errors.
		prec::estimate("funcName", f_approx, f_exact, options);

	// Stop testing the module and write the results
	prec::terminate();
}
