
///
/// @file statistics.cpp Basic statistical functions on a single data set.
/// This example may be compiled using 'make statistics'
///

#include "theoretica.h"

// The utility.h header is not included
// by theoretica.h and is needed for insert_data()
#include "utility.h"
using namespace th;


int main(int argc, char const *argv[]) {


	std::vector<real> X;

	// Get the data from standard input
	std::cout << "Insert X (write END to stop):" << std::endl;
	insert_data(X, "END");


	// You can easily compute statistical functions on the data set
	// by using functionalities defined in "statistics.h"
	std::cout << std::endl;
	std::cout << "N = " << X.size() << std::endl;
	std::cout << "Mean: " << mean(X) << std::endl;
	std::cout << "Variance: " << sample_variance(X) << std::endl;

	// Many functions with long names have short-hands:
	// sample_standard_deviation = smpl_stdev
	// sample_mean_standard_deviation = smpl_stdom
	std::cout << "Standard Deviation: " << smpl_stdev(X) << std::endl;
	std::cout << "Mean Standard Deviation: " << smpl_stdom(X) << std::endl;
	
	std::cout << "Relative Error: " << sample_standard_relative_error(X) * 100 << "%" << std::endl;

	return 0;
}
