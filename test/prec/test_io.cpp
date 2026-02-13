
///
/// @file test_io.cpp Input and output test unit.
///

#include "theoretica.h"
#include "chebyshev.h"

#include "io/format_csv.h"

using namespace chebyshev;
using namespace theoretica;


template<unsigned int N>
real absmax(const vec<real, N>& v) {

	real max = 0.0;
	for (real x : v)
		max = std::max(max, std::abs(x));

	return max;
}


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("io");
	ctx.settings.outputFiles = { "test/prec/prec_io.csv" };
	random::random_source rnd = ctx.random->get_rnd();
	
	io::println("If you see this, everything is going as intended.");

	// Vector to CSV without header
	{
		unsigned int N = 1E+03;

		// Generate a random vector
		vec<real> v (N);
		for (size_t i = 0; i < N; i++)
			v[i] = rnd.gaussian(0, 1);
		
		// Back-and-forth write and reading test
		io::write_csv("./test/prec/test.csv", v);

		vec<real> w;
		io::read_csv("./test/prec/test.csv", w);

		ctx.equals("write_csv/read_csv(vec<real>)", absmax(v - w), 0.0, 1E-07);
	}

	// Vector to CSV without header
	{

		// Generate a random vector
		vec<real, 1000> v;
		for (size_t i = 0; i < v.size(); i++)
			v[i] = rnd.gaussian(0, 1);
		
		// Back-and-forth write and reading test
		io::write_csv("./test/prec/test.csv", v);

		vec<real, 1000> w;
		io::read_csv("./test/prec/test.csv", w);

		ctx.equals("write_csv/read_csv(vec<real, N>)", absmax(v - w), 0.0, 1E-07);


		// Check that, for mismatched sizes, the remaining empty elements are filled with NaN
		vec<real, 1001> z;
		io::read_csv("./test/prec/test.csv", z);

		ctx.equals("read_csv(vec<real, N>) (NaN)", is_nan(z[1000]), true, 0);
	}


	// Vector to CSV with header
	{
		unsigned int N = 1E+03;

		// Generate a random vector
		vec<real> v (N);
		for (size_t i = 0; i < N; i++)
			v[i] = rnd.gaussian(0, 1);
		
		// Back-and-forth write and reading test
		io::write_csv("./test/prec/test.csv", v, "Vector");

		vec<real> w;
		io::read_csv("./test/prec/test.csv", w);

		ctx.equals("write_csv/read_csv(vec<real>, header)", absmax(v - w), 0.0, 1E-07);
	}
}
