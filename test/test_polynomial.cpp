
///
/// @file test_polynomial.cpp Polynomial class test cases
///

#include "theoretica.h"
#include <cmath>
#include "chebyshev/prec.h"

using namespace chebyshev;
using namespace theoretica;


prec::estimate_result test_polynomial_eval(interval k, Real tol, unsigned int n) {
	
	Real max = 0, sum = 0, sum2 = 0;

	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (unsigned int i = 0; i < n; ++i) {
		
		unsigned int degree = g() % 100;
		polynomial<real> P = polynomial<>();
		P.coeff.reserve(degree + 1);

		real s = 0;

		for (unsigned int j = 0; j < degree + 1; ++j) {
			P.coeff.push_back(rand_uniform(k.a, k.b, g));
			s += P[j];
		}

		real diff = th::abs(P(1.0) - s);

		if(max < diff)
			max = diff;

		sum += diff;
		sum2 += square(diff);
	}

	prec::estimate_result p;
	p.max_err = max;
	p.abs_err = sum / n;
	p.rms_err = th::sqrt(sum2) / n;
	p.mean_err = sum / n;

	// Undefined relative error
	p.rel_err = 0;

	if(p.max_err > tol)
		p.failed = true;

	return p;
}


prec::estimate_result test_polynomial_div(interval k, Real tol, unsigned int n) {
	
	real max = 0, sum = 0, sum2 = 0;
	prec::estimate_result p;
	PRNG g = PRNG::xoshiro(time(nullptr));
	g.discard(1000);

	for (unsigned int i = 0; i < n; ++i) {

		unsigned int A_degree = 1 + (g() % 5);
		unsigned int B_degree = 1 + (g() % 5);

		polynomial<real> A = polynomial<>(), B = polynomial<>(), C = polynomial<>();

		for (unsigned int j = 0; j < A_degree + 1; ++j)
			A.coeff.push_back(rand_uniform(k.a, k.b, g));

		for (unsigned int j = 0; j < B_degree + 1; ++j)
			B.coeff.push_back(rand_uniform(k.a, k.b, g));

		C = A * B;
		polynomial<real> res = C / B;
		const unsigned int res_degree = res.find_order();

		if(A_degree != res_degree) {
			p.failed = true;
			std::cout << "\t\tFAILED: C = A*B = " << C << std::endl;
		}

		for (unsigned int j = 0; j < min(A.size(), res.size()); ++j) {
			
			real diff = th::abs(A[j] - res[j]);
			sum += diff / res_degree;
			sum2 += square(diff / res_degree);

			if(max < diff)
				max = diff;
		}	
	}

	p.max_err = max;
	p.abs_err = sum / n;
	p.rms_err = th::sqrt(sum2) / n;
	p.mean_err = sum / n;

	// Undefined relative error
	p.rel_err = 0;

	if(p.max_err > tol)
		p.failed = true;

	return p;
}


int main(int argc, char const *argv[]) {

	const real MAX = 1000000;
	const real MIN = -MAX;

	prec::state.outputFolder = "test/";
	
	prec::setup("polynomial");

		prec::estimate("polynomial<>::eval", test_polynomial_eval, interval(MIN, MAX));


		// Investigate polynomial division algorithm not converging

		// prec::estimate("polynomial<>::operator/", test_polynomial_div, interval(-100, 100),
		// 	prec::state.defaultTolerance, false, 20);

	prec::terminate();
}
