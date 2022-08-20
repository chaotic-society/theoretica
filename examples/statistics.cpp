
///
/// @file statistics.cpp Basic statistical functions on a single data set
///

#include "theoretica.h"

// The utility.h header is not included
// by theoretica.h and is needed for insert_data()
#include "utility.h"
using namespace th;


int main(int argc, char const *argv[]) {

	std::cout.precision(8);

	vec_buff X;

	std::cout << "Insert X (write END to stop):" << std::endl;
	insert_data(X, "END");

	std::cout << std::endl;
	std::cout << "N = " << X.size() << std::endl;
	std::cout << "Mean: " << mean(X) << std::endl;
	std::cout << "Variance: " << sample_variance(X) << std::endl;
	std::cout << "Standard Deviation: " << smpl_stdev(X) << std::endl;
	std::cout << "Relative Error: " << sample_standard_relative_error(X) * 100 << "%" << std::endl;
	std::cout << "Mean Standard Deviation: " << smpl_stdom(X) << std::endl;

	return 0;
}
