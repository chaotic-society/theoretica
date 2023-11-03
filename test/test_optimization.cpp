
#include "theoretica.h"
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


// Multivariate root search
// {1, e} is a root
template<typename NumType>
inline vec<NumType, 2> f1(vec<NumType, 2> v) {

	const NumType x = v[0];
	const NumType y = v[1];

	return {
		exp(x) - y,
		x * y - exp(x)
	};
}


int main(int argc, char const *argv[]) {

	prec::state.outputFolder = "test/";
	
	prec::setup("optimization");


		prec::equals(
			"multiroot_newton (1)",
			multiroot_newton<2>(f1, {2, 2}, 10E-10)[0],
			1
		);


		prec::equals(
			"multiroot_newton (2)",
			multiroot_newton<2>(f1, {2, 2}, 10E-10)[1],
			E
		);


	prec::terminate();
}
