
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	prec::setup("signal", argc, argv);

		output::settings.outputFiles = { "test/prec/prec_signal.csv" };


		// Test fft.h

	{
		cvec x = {};
		cvec empty = {};

		prec::equals(
			"fft({})",
			signal::fft(x) == empty,
			true, 0
		);
	}

	{
		cvec x = {};
		cvec empty = {};

		prec::equals(
			"ifft({})",
			signal::ifft(x) == empty,
			true, 0
		);
	}

	{
		PRNG g = PRNG::xoshiro(time(nullptr));
		pdf_sampler gauss = pdf_sampler::gaussian(0, 1E+03, g);


		unsigned int N = (1 << 16);
		cvec x = cvec(N);
		gauss.fill(x);

		prec::equals(
			"ifft(fft(x)) = x",
			algebra::linf_norm(signal::ifft(signal::fft(x)) - x),
			0
		);
	}

	{
		cvec x = {1, 1};
		cvec expected = {2, 0};

		prec::equals(
			"fft(1, 1)",
			algebra::linf_norm(signal::fft(x) - expected),
			0
		);
	}

	{
		cvec x = {1, 1, 1};

		prec::equals(
			"fft (N != 2^m)",
			is_nan(signal::fft(x)[0]),
			true
		);
	}

	prec::terminate();
}
