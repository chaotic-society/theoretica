
#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


int main(int argc, char const *argv[]) {
	
	prec::setup("signal");

		output::state.outputFiles = { "test/prec_signal.csv" };

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

	prec::terminate();
}
