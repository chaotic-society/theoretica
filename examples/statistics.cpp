
///
/// @file statistics.cpp Basic statistical functions on a single data set.
/// This example may be compiled using 'make statistics'
///

#include "theoretica.h"

// The utility.h header is not included
// by theoretica.h and is needed for insert_data()
#include "utility.h"
using namespace th;


int main() {

	vec<real> X;

	// Get the data from standard input
	println("Insert dataset (empty line to stop):");
	X = readln();


	// You can easily compute statistical functions on the data set
	// by using functions defined in "statistics.h"
	println();
	println("N =", X.size());
	println("Mean:", stats::mean(X));
	println("Variance:", stats::variance(X));
	println("Standard Deviation:", stats::stdev(X));
	println("Mean Standard Deviation:", stats::stdom(X));
	println("Relative Error:", stats::standard_relative_error(X) * 100, "%");
}
