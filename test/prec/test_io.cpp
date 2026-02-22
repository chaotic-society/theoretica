
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

	auto vec_opt = prec::equation_options<vec<real>>(
		1E-08, prec::distance::euclidean<vec<real>>
	);

	{
		ctx.equals("is_number", io::is_number("3.1415"), true);
		ctx.equals("is_number", io::is_number("1,414"), true);
		ctx.equals("is_number", io::is_number("123"), true);
		ctx.equals("is_number", io::is_number("NaN"), true);
		ctx.equals("is_number", io::is_number("nan"), true);
		ctx.equals("is_number", io::is_number("+inf"), true);
		ctx.equals("is_number", io::is_number("-1.0E+99"), true);
		ctx.equals("is_number", io::is_number(""), false);
		ctx.equals("is_number", io::is_number("Hello, World!"), false);
		ctx.equals("is_number", io::is_number("Clearly not a number, but not a Not a Number"), false);
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
		io::write_csv("./test/prec/test.csv", "Vector", v);

		vec<real> w;
		io::read_csv("./test/prec/test.csv", w);

		vec<real> z;
		io::read_csv("./test/prec/test.csv", "Vector", z);

		ctx.equals("write_csv/read_csv(vec<real>, header)", absmax(v - w), 0.0, 1E-07);
		ctx.equals("read_csv(vec<real>, col_name)", absmax(v - z), 0.0, 1E-07);
	}


	// Matrix to CSV
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

	// Histogram to CSV
	{
		vec<real> v (1000);
		rnd.gaussian(v, 0, 1);

		histogram hist (v);
		io::write_csv("./test/prec/test.csv", hist);
		
		histogram hist2;
		io::read_csv("./test/prec/test.csv", hist2);

		auto bin_opt = prec::equation_options<std::vector<unsigned int>>(
			1E-08, prec::distance::euclidean<std::vector<unsigned int>>
		);
		
		ctx.equals("write_csv/read_csv(histogram)", hist.bins(), hist2.bins(), bin_opt);
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

	// Table creation and insertion
	{
		data_table table;

		ctx.equals("data_table.empty()", table.empty(), true);

		table.insert("A", {1, 2, 3});
		table.insert("B", {4, 5});
		table.insert("C", {6, 7, 8, 9});

		ctx.equals("data_table.columns()", table.cols(), 3, 0);
		ctx.equals("data_table.rows()", table.rows(), 4, 0);
		ctx.equals("data_table[\"A\"][1]", table["A"][1], 2.0, 0);
	}

	// Table creation from a hashmap
	{
		vec<real> v = {1, 2, 3};
		vec<real> w = {th::PI, th::E};

		std::map<std::string, vec<real>> m;
		m["v"] = v;
		m["w"] = w;

		data_table table1 (m);
		ctx.equals("data_table(map)", table1["v"], v, vec_opt);
		ctx.equals("data_table(map)", table1["w"], w, vec_opt);

		data_table table2 = data_table(3, {"v1", "v2"});
		ctx.equals("data_table(n_rows, col_name)", table2["v1"].size(), 3);
		ctx.equals("data_table(n_rows, col_name)", table2["v2"].size(), 3);
	}
}
