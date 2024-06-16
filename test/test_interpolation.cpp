
#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;



int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("interpolation");

		

	prec::terminate();
}
