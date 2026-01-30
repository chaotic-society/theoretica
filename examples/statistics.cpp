
///
/// @file statistics.cpp Basic statistical functions on a single data set.
/// This example may be compiled using 'make statistics'
///

#include "theoretica.h"
#include "io/io.h"
using namespace th;


int main() {

	vec<real> X;

	// Get the data from standard input
	io::println("Insert dataset (empty line to stop):");

	std::string line = "a";
	while (line != "") {
		
		line = io::readln();

		if(line != "") {
			real value = std::stod(line);
			X.append(value);
		}
	}


	// You can easily compute statistical functions on the data set
	// by using functions defined in "statistics.h"
	io::println();
	io::println("N =", X.size());
	io::println("Mean:", stats::mean(X));
	io::println("Variance:", stats::variance(X));
	io::println("Standard Deviation:", stats::stdev(X));
	io::println("Mean Standard Deviation:", stats::stdom(X));
	io::println("Relative Error:", stats::standard_relative_error(X) * 100, "%");
}
