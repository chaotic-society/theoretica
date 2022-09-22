
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


prec::estimate_result test_custom(interval k, Real tol, unsigned int n) {
	// Custom test
}


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("module name");

		prec::estimate(...);

		prec::equals(...);

	prec::terminate();
}
