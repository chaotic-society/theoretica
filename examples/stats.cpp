
///
/// @file stats.cpp Basic statistical functions on a data set
///

#include "uroboro.h"

// The utility.h header is not included
// by uroboro.h
#include "utility.h"
using namespace umath;


int main(int argc, char const *argv[]) {

	vec_buff data;

	std::cout << "Insert data (write END to stop):" << std::endl;
	insert_data(data, "");

	std::cout << std::endl;
	print_sample_stats(data);

	return 0;
}
