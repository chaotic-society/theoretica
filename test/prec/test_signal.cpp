
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("signal", argc, argv);
	ctx.output->settings.outputFiles = { "test/prec/prec_signal.csv" };

	// Test fft.h

	// On empty vector, return { NaN }
	{
		cvec x;
		cvec result;

		try {
			result = signal::fft(x);
		} catch (std::exception& e) {
			result = {nan()};
		}

		ctx.equals(
			"fft({})",
			result.size() == 1 && is_nan(result[0]),
			true, 0
		);
	}

	{
		cvec x;
		cvec result;

		try {
			result = signal::ifft(x);
		} catch (std::exception& e) {
			result = {nan()};
		}

		ctx.equals(
			"ifft({})",
			result.size() == 1 && is_nan(result[0]),
			true, 0
		);
	}

	{
		PRNG g = PRNG::xoshiro(time(nullptr));
		pdf_sampler gauss = pdf_sampler::gaussian(0, 1E+03, g);


		unsigned int N = (1 << 16);
		cvec x = cvec(N);
		gauss.fill(x);

		ctx.equals(
			"ifft(fft(x)) = x",
			algebra::linf_norm(signal::ifft(signal::fft(x)) - x),
			0
		);
	}

	{
		cvec x = {1, 1};
		cvec expected = {2, 0};

		ctx.equals(
			"fft(1, 1)",
			algebra::linf_norm(signal::fft(x) - expected),
			0
		);
	}

	{
		cvec x = {1, 1, 1};

		ctx.equals(
			"fft (N != 2^m)",
			is_nan(signal::fft(x)[0]),
			true
		);
	}
}
