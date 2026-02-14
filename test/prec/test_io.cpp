
///
/// @file test_io.cpp Input and output test unit.
///

#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


template<unsigned int N>
real absmax(const vec<real, N>& v) {

	real max = 0.0;
	for (real x : v) {

		const real abs_x = std::abs(x);
		
		// Check for NaN
		if (abs_x != abs_x)
			max = inf();

		max = std::max(max, abs_x);
	}

	return max;
}

template<unsigned int N, unsigned int K>
real absmax(const mat<real, N, K>& v) {

	real max = 0.0;
	for (real x : v) {

		const real abs_x = std::abs(x);
		
		if (abs_x != abs_x)
			max = inf();

		max = std::max(max, abs_x);
	}

	return max;
}


int main(int argc, char const *argv[]) {
	
	auto ctx = prec::make_context("io");
	ctx.settings.outputFiles = { "test/prec/prec_io.csv" };
	random::random_source rnd = ctx.random->get_rnd();

	// io.h
	
	io::println("If you see this, everything is going as intended.");


	// strings.h

	auto str_opt = prec::equation_options<std::string>(
		0.0, prec::distance::hamming
	);

	{
		ctx.equals("trim", io::trim("  \"Hello, World!\"  "), std::string("\"Hello, World!\""), str_opt);
		ctx.equals("trim", io::trim(" 	   many words here  !!!   	 "), std::string("many words here  !!!"), str_opt);
		ctx.equals("unquote", io::unquote("\"Hello, World!\""), std::string("Hello, World!"), str_opt);
		ctx.equals("unquote", io::unquote("\"this quote is open"), std::string("\"this quote is open"), str_opt);
	}


	// format_csv.h

	// Vector to CSV without header
	{
		unsigned int N = 1E+03;

		// Generate a random vector
		vec<real> v (N);
		rnd.gaussian(v, 0, 1);
		
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
		rnd.gaussian(v, 0, 1);
		
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
		rnd.gaussian(v, 0, 1);
		
		// Back-and-forth write and reading test
		io::write_csv("./test/prec/test.csv", v, "Vector");

		vec<real> w;
		io::read_csv("./test/prec/test.csv", w);

		ctx.equals("write_csv/read_csv(vec<real>, header)", absmax(v - w), 0.0, 1E-07);
	}


	// Matrix to CSV without header
	{

		// Generate a random matrix
		mat<real, 100, 100> A;
		rnd.gaussian(A, 0, 1);

		// Back-and-forth write and reading test
		io::write_csv("./test/prec/test.csv", A);

		mat<real> B;
		io::read_csv("./test/prec/test.csv", B);
		ctx.equals("write_csv/read_csv(mat<real, N, M>)", absmax(A - B), 0.0, 1E-07);

		mat<real> C = io::read_csv<mat<real>>("./test/prec/test.csv");
		ctx.equals("write_csv/read_csv<mat<real>>()", absmax(A - C), 0.0, 1E-07);
	}

	// Test CSV tokenization and parsing
	{

		std::string line = "  1.2,  \"3151,	 726\", 	  \"135.153161,135136\"   	";
		std::vector<std::string> tokens = io::parse_csv(line);

		ctx.equals("parse_csv", tokens.size(), 3, 0);

		if (tokens.size() >= 3) {
			ctx.equals("parse_csv", tokens[0], std::string("1.2"), str_opt);
			ctx.equals("parse_csv", tokens[1], std::string("3151,	 726"), str_opt);
			ctx.equals("parse_csv", tokens[2], std::string("135.153161,135136"), str_opt);
		}
	}


	// data_table.h

	{
		data_table table;

		table.insert("A", {1, 2, 3});
		table.insert("B", {4, 5});
		table.insert("C", {6, 7, 8, 9});

		ctx.equals("data_table.columns()", table.cols(), 3, 0);
		ctx.equals("data_table.rows()", table.rows(), 4, 0);
		ctx.equals("data_table[\"A\"][1]", table["A"][1], 2.0, 0);
	}
}
