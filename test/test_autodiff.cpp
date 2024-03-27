
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("autodiff");

		// Compare the automatic derivative to the analytical derivative

		prec::estimate("dual::Dual()",
			[](real x) {
				dual d = dual(x, 1);
				return (th::cos(square(d)) / th::exp(-square(d)) / ln(1 / square(d))).Dual();
			},
			[](real x) {
				return (2 * th::exp(square(x)) * ((square(x) * ln(1 / square(x)) + 1)
							* th::cos(square(x)) - square(x) * ln(1 / square(x)) * th::sin(square(x))))
								/ (x * square(ln(1 / square(x))));
			}, {interval(0.001, 0.5), interval(-0.5, -0.001)});	

	prec::terminate();
}
