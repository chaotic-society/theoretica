
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("module name");

		prec::estimate(...);

		prec::equals(...);

	prec::terminate();
}
