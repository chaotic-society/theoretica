
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
	// by using functionalities defined in "statistics.h"
	println();
	println("N =", X.size());
	println("Mean:", mean(X));
	println("Variance:", sample_variance(X));

	// Many functions with long names have short-hands:
	// sample_standard_deviation = smpl_stdev
	// sample_mean_standard_deviation = smpl_stdom
	println("Standard Deviation:", smpl_stdev(X));
	println("Mean Standard Deviation:", smpl_stdom(X));
	
	println("Relative Error:", sample_standard_relative_error(X) * 100, "%");

	return 0;
}
