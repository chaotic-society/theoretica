#define UROBORO_INCLUDE_ALL
#include "./statistics.h"
#include "./utility.h"
#include <iostream>
#include <string>

using namespace uroboro;


int main(int argc, char const *argv[]) {

	std::cout << "How many records to insert? (1 or 2 only)" << std::endl;
	std::string line;
	std::getline(std::cin, line);

	bool two_records;

	if(line == "1")
		two_records = false;
	else if(line == "2")
		two_records	= true;
	else {
		std::cout << "Insertion error" << std::endl;
		exit(1);
	}

	std::cout << "\nInsert each value and press Enter (END to stop insertion)" << std::endl;

	vec_buff X;
	vec_buff Y;

	std::cout << "Insert X values:" << std::endl;
	insert_data(X, "");

	if(two_records) {
		std::cout << "Insert Y values:" << std::endl;
		insert_data(Y, "");
	}

	std::cout << "\nStats for X:" << std::endl;
	print_sample_stats(X);

	if(two_records) {
		std::cout << "\nStats for Y:" << std::endl;
		print_sample_stats(Y);
	}

	if(two_records) {

		std::cout << "Covariance: " << sample_covariance(X, Y) << "\n\n";

		real intercept = least_squares_linear_intercept(X, Y);
		real slope = least_squares_linear_slope(X, Y);

		std::cout << "Minimum Squares Linearization:\n";
		std::cout << "Minimum Squares Intercept: " << intercept << std::endl;
		std::cout << "Minimum Squares Slope: " << slope << std::endl;
		std::cout << "Minimum Squares Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << std::endl;

		std::cout << "Correlation Coefficient: " <<
			sample_correlation_coefficient(X, Y) << std::endl;
		std::cout << std::endl;
	}

	std::cout << "Press Enter to exit..." << std::endl;
	std::cin.get();

	return 0;
}
