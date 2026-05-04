
///
/// @file test_io.cpp Input and output test unit.
///

#include "theoretica.h"
#include "chebyshev.h"

#ifdef THEORETICA_HAS_HDF5
#include "io/hdf5.h"
#endif

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


	// csv.h

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
		ctx.equals("data_table.size()", table.size(), 0, 0);

		table.insert("A", {1, 2, 3});
		table.insert("B", {4, 5});
		table.insert("C", {6, 7, 8, 9});

		ctx.equals("data_table.header()", table.header() == std::vector<std::string>({"A", "B", "C"}), true);
		ctx.equals("data_table.has_column()", table.has_column("A"), true);
		ctx.equals("data_table.has_column()", table.has_column("D"), false);

		ctx.equals("data_table.cols()", table.cols(), 3, 0);
		ctx.equals("data_table.rows()", table.rows(), 4, 0);
		ctx.equals("data_table.size()", table.size(), 9, 0);
		ctx.equals("data_table[\"column\"][1]", table["A"][1], 2.0, 0);

		ctx.equals("data_table.at(string, idx)", table.at("B", 0), 4.0, 0);
		ctx.equals("data_table.row(idx)", table.row(2)["C"], 8.0, 0);
		ctx.equals("data_table.row(idx)", is_nan(table.row(2)["B"]), true, 0);
	}

	// Table creation from hashmap
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

	// Table manipulation
	{
		data_table table = data_table(4, {"length", "value", "count"});
		table["length"] = {1, 2, 3, 4};
		table["value"] = {11, 22, 31, 45};
		table["count"] = {1, 3, 5, 3};

		table.rename("count", "number");
		ctx.equals("data_table.rename()", table["number"], vec<real>({1, 3, 5, 3}), vec_opt);
		ctx.equals("data_table.rename()", table.has_column("count"), false);

		auto new_table = table.select({"length", "value"});
		ctx.equals("data_table.select()", new_table.header() == std::vector<std::string>({"length", "value"}), true);
		ctx.equals("data_table.select()", new_table["length"], vec<real>({1, 2, 3, 4}), vec_opt);

		table.drop_column("v");
		ctx.equals("data_table.drop_column()", table.has_column("v"), false);

		table.insert("v", vec<real>({7, 8, 9}));
		table.drop_columns({"v", "w"});
		ctx.equals("data_table.drop_columns()", table.has_column("v"), false);
		ctx.equals("data_table.drop_columns()", table.has_column("w"), false);

		table.clear();
		ctx.equals("data_table.clear()", table.empty(), true);
	}

	// Table data handling
	{
		data_table table;
		vec<real> v = vec<real>({3, 4, 5});
		table.insert("v", v);

		ctx.equals("data_table.data()", table.data().size(), 1);
		ctx.equals("data_table.data()", table.data()[0], v, vec_opt);
		ctx.equals("data_table[0]", table[0], v, vec_opt);

		table.insert("w", vec<real>({1, 2, 3}));

		auto table_head = table.head(2);
		ctx.equals("data_table.head()", table_head.rows(), 2);
		ctx.equals("data_table.head()", table_head["v"], vec<real>({3, 4}), vec_opt);
		ctx.equals("data_table.head()", table_head["w"], vec<real>({1, 2}), vec_opt);

		auto table_tail = table.tail(1);
		ctx.equals("data_table.tail()", table_tail.rows(), 1);
		ctx.equals("data_table.tail()", table_tail["v"], vec<real>({5}), vec_opt);
		ctx.equals("data_table.tail()", table_tail["w"], vec<real>({3}), vec_opt);

		mat<real> A = table.to_matrix();
		mat<real> A_expected = {{3, 1}, {4, 2}, {5, 3}};
		ctx.equals("data_table.to_matrix()", absmax(A - A_expected), 0);

		mat<real> B = {{1, -1}, {2, -2}, {3, -3}};
		table.from_matrix(B, {"v", "w"});
		ctx.equals("data_table.from_matrix()", table["v"], vec<real>({1, 2, 3}), vec_opt);
		ctx.equals("data_table.from_matrix()", table["w"], vec<real>({-1, -2, -3}), vec_opt);
	}


	// hdf5.h

#ifdef THEORETICA_HAS_HDF5

	vec<real> v (1000);
	rnd.gaussian(v, 0, 1E+09);

	mat<real> A (1000, 1000);
	rnd.gaussian(A, 0, 1E+09);

	io::println("Writing to HDF5 file...");
	// Write to file
	{	
		io::hdf5_file file ("./test/prec/test.h5", true);

		try {
			file.create_group("/group");
		} catch (const std::exception& e) {
			ctx.equals("hdf5_file.create_group()", false, true);
		}

		try {
			file.write_vec("/group/vec", v);
		} catch (const std::exception& e) {
			ctx.equals("hdf5_file.write_vec()", false, true);
		}

		try {
			file.write_mat("/group/mat", A);
		} catch (const std::exception& e) {
			ctx.equals("hdf5_file.write_mat()", false, true);
		}

		try {
			file.write_attribute("/group/vec", "author", "Albert Einstein");
		} catch (const std::exception& e) {
			ctx.equals("hdf5_file.write_attribute()", false, true);
		}
	}

	io::println("Reading from HDF5 file...");
	// Read from file
	{
		io::hdf5_file file ("./test/prec/test.h5", false);

		vec<real> v_read = file.read_vec("/group/vec");
		ctx.equals("hdf5_file.read_vec()", absmax(v - v_read), 0.0);

		mat<real> A_read = file.read_mat("/group/mat");
		ctx.equals("hdf5_file.read_mat()", absmax(A - A_read), 0.0);

		std::string author = file.read_attribute<std::string>("/group/vec", "author");
		ctx.equals("hdf5_file.read_attribute()", author, std::string("Albert Einstein"), str_opt);
	
		ctx.equals("hdf5_node.is_group()", file["group"].is_group(), true);
		ctx.equals("hdf5_node.is_dataset()", file["group"]["vec"].is_dataset(), true);

	}

#endif

}
