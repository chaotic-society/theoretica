#define UROBORO_INCLUDE_ALL
#include "./statistics.h"
#include "./utility.h"
#include <iostream>
#include <string>

using namespace uroboro;


int main(int argc, char const *argv[]) {

	std::cout << "1:\tInsert 1 record" << std::endl;
	std::cout << "2:\tInsert 2 records" << std::endl;
	std::cout << "3:\tInsert 3 records" << std::endl;
	std::cout << "4:\tPropagate sum" << std::endl;
	std::cout << "5:\tPropagate product (or quotient)" << std::endl;

	std::string line;
	size_t n_records;

	std::getline(std::cin, line);

	if(line == "1")
		n_records = 1;
	else if(line == "2")
		n_records	= 2;
	else if(line == "3")
		n_records	= 3;
	else if(line == "4")
		n_records	= 4;
	else if(line == "5")
		n_records	= 5;
	else {
		std::cout << "Input error\n";
		exit(1);
	}

	if(n_records > 5) {
		std::cout << "Input error" << std::endl;
		exit(1);
	}

	std::cout << "\nInsert each value and press Enter (END to stop insertion)" << std::endl;

	vec_buff X, Y, W;

	std::cout << "Insert X values:" << std::endl;
	insert_data(X, "");

	if(n_records == 4) {
		std::cout << "Propagated error: " << propagate_sum(X) << std::endl;
		std::cout << "(X values used as stdev of variables)" << std::endl;
		exit(0);
	}

	if(n_records > 1) {
		std::cout << "Insert Y values:" << std::endl;
		insert_data(Y, "");
	}

	if(n_records == 5) {
		std::cout << "Propagated error: " << propagate_product(X, Y) << std::endl;
		std::cout << "(X values used as stdev of variables)" << std::endl;
		std::cout << "(Y values used as mean of variables)" << std::endl;
		exit(0);
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

		if(X.size() != Y.size()) {
			std::cout << "ERROR: Data sets must have the same size" << std::endl;
			exit(1);
		}

		vec_buff weights = Y;
		for (int i = 0; i < weights.size(); ++i) {
			weights[i] = 1.0 / square(weights[i]);
		}

		real r = sample_correlation_coefficient(X, Y);

		std::cout << "Covariance: " << sample_covariance(X, Y) << "\n";
		std::cout << "Correlation Coefficient: " << r << "\n";
		std::cout << "r-Squared: " << square(r) << "\n\n";

		real intercept = least_squares_linear_intercept(X, Y);
		real slope = least_squares_linear_slope(X, Y);

		std::cout << "MINIMUM SQUARES LINEARIZATION:\n";
		std::cout << "y = A + Bx\n";
		std::cout << "Minimum Squares Intercept: " << intercept << "\n";
		std::cout << "Minimum Squares Slope: " << slope << "\n";
		std::cout << "Minimum Squares Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << "\n";
		std::cout << std::endl;

		std::cout << "Weighted mean: " << weighted_mean(X, weights) << std::endl;
		std::cout << "Weighted mean sigma: " << (1.0 / sqrt(sum(weights))) << std::endl;
		std::cout << std::endl;
	}


	if(n_records == 3) {

		if(X.size() != W.size()) {
			std::cout << "ERROR: Data sets must have the same size" << std::endl;
			exit(1);
		}

		vec_buff sigma = W;

		// Calculate weights as 1 / sigma^2
		for (int i = 0; i < W.size(); ++i) {
			W[i] = 1.0 / square(sigma[i]);
		}

		real intercept = least_squares_weighted_linear_intercept(X, Y, W);
		real slope = least_squares_weighted_linear_slope(X, Y, W);

		std::cout << "WEIGHTED MINIMUM SQUARES LINEARIZATION:\n";
		std::cout << "Weighted Minimum Squares Intercept: " << intercept << std::endl;
		std::cout << "Weighted Minimum Squares Slope: " << slope << std::endl;
		std::cout << "Weighted Minimum Squares Error: " <<
			least_squares_linear_error(X, Y, intercept, slope) << std::endl;
		std::cout << std::endl;

		std::cout << "Minimum Squares Sigma A (sigma_y = w[0]): " <<
			least_squares_linear_sigma_A(X, Y, sigma[0]) << std::endl;

		std::cout << "Minimum Squares Sigma B (sigma_y = w[0]): " <<
			least_squares_linear_sigma_B(X, Y, sigma[0]) << "\n" << std::endl;

		std::cout << "CHI-SQUARE on WEIGHTED LINEARIZATION: " <<
			chi_square_linearization(X, Y, sigma, intercept, slope) << std::endl;

		std::cout << "REDUCED CHI-SQUARE on WEIGHTED LINEARIZATION: " <<
			reduced_chi_square_linearization(X, Y, sigma, intercept, slope) << std::endl;
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

		std::cout << "Log Covariance: " << sample_covariance(X, Y) << std::endl;

		std::cout << "Log Correlation Coefficient: " <<
			sample_correlation_coefficient(X, Y) << std::endl;

	}

	std::cout << "Press Enter to exit..." << std::endl;
	std::cin.get();

	return 0;
}
