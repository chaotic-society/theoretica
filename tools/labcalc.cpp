#define UROBORO_INCLUDE_ALL
#include "./statistics.h"
#include "./utility.h"
#include <iostream>
#include <string>

using namespace uroboro;


int main(int argc, char const *argv[]) {

	std::cout << "How many records to insert? (1, 2 or 3)" << std::endl;
	std::string line;
	std::getline(std::cin, line);

	size_t n_records;

	if(line == "1")
		n_records = 1;
	else if(line == "2")
		n_records	= 2;
	else if(line == "3")
		n_records	= 3;
	else {
		std::cout << "Insertion error\n";
		exit(1);
	}

	std::cout << "\nInsert each value and press Enter (END to stop insertion)\n";

	vec_buff X, Y, W;

	std::cout << "Insert X values:" << std::endl;
	insert_data(X, "");

	if(n_records > 1) {
		std::cout << "Insert Y values:" << std::endl;
		insert_data(Y, "");
	}

	if(n_records > 2) {
		std::cout << "Insert W values:" << std::endl;
		insert_data(W, "");
	}

	if(n_records == 1) {
		std::cout << "\nStats for X:" << std::endl;
		print_sample_stats(X);
	}


	if(n_records > 1) {

		if(X.size() != Y.size() || X.size() != W.size()) {
			std::cout << "ERROR: Data sets must have the same size" << std::endl;
			exit(1);
		}

		std::cout << "Covariance: " << sample_covariance(X, Y) << "\n";
		std::cout << "Correlation Coefficient: " <<
			sample_correlation_coefficient(X, Y) << "\n\n";

		real intercept = least_squares_linear_intercept(X, Y);
		real slope = least_squares_linear_slope(X, Y);

		std::cout << "MINIMUM SQUARES LINEARIZATION:\n";
		std::cout << "y = A + Bx\n";
		std::cout << "Minimum Squares Intercept: " << intercept << "\n";
		std::cout << "Minimum Squares Slope: " << slope << "\n";
		std::cout << "Minimum Squares Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << "\n";
		std::cout << std::endl;
	}


	if(n_records == 3) {

		// Calculate weights as 1 / sigma^2
		for (int i = 0; i < W.size(); ++i) {
			W[i] = 1.0 / square(W[i]);
		}

		real intercept = least_squares_weighted_linear_intercept(X, Y, W);
		real slope = least_squares_weighted_linear_slope(X, Y, W);

		std::cout << "WEIGHTED MINIMUM SQUARES LINEARIZATION:\n";
		std::cout << "Weighted Minimum Squares Intercept: " << intercept << std::endl;
		std::cout << "Weighted Minimum Squares Slope: " << slope << std::endl;
		std::cout << "Weighted Minimum Squares Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << std::endl;
		std::cout << std::endl;

		for (int i = 0; i < X.size(); ++i) {
			X[i] = uroboro::ln(X[i]);
		}

		for (int i = 0; i < Y.size(); ++i) {
			Y[i] = uroboro::ln(Y[i]);
		}

		intercept = least_squares_weighted_linear_intercept(X, Y, W);
		slope = least_squares_weighted_linear_slope(X, Y, W);

		std::cout << "WEIGHTED MINIMUM SQUARES LOGARITHM LINEARIZATION:\n";
		std::cout << "ln(y) = A + Bln(x)" << std::endl;
		std::cout << "Weighted Minimum Squares Log Intercept: " << intercept << std::endl;
		std::cout << "Weighted Minimum Squares Log Slope: " << slope << std::endl;
		std::cout << "Weighted Minimum Squares Log Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << std::endl;
		std::cout << std::endl;

	}

	std::cout << "Press Enter to exit..." << std::endl;
	std::cin.get();

	return 0;
}
