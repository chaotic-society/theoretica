
///
/// @file stats.cpp Basic statistical functions on a data set
///

#include "theoretica.h"

// The utility.h header is not included
// by theoretica.h
#include "utility.h"
using namespace th;


int main(int argc, char const *argv[]) {

	std::cout.precision(8);

	int n;
	vec_buff X;
	vec_buff Y;
	vec_buff Z;
	real sigma_Y;

	std::cout << "1: Single dataset\n" << "2: Two datasets\n" << "3: Three datasets\n";
	std::cin >> n;

	if(n == 1) {

		std::cout << "Insert X (write END to stop):" << std::endl;
		insert_data(X, "END");
		std::cout << std::endl;
		std::cout << "N = " << X.size() << std::endl;
		std::cout << "Mean: " << mean(X) << std::endl;
		std::cout << "Variance: " << sample_variance(X) << std::endl;
		std::cout << "Standard Deviation: " << smpl_stdev(X) << std::endl;
		std::cout << "Relative Error: " << sample_standard_relative_error(X) * 100 << "%" << std::endl;
		std::cout << "Mean Standard Deviation: " << smpl_stdom(X) << std::endl;

	} else if(n == 2) {

		std::cout << "Insert X (write END to stop):" << std::endl;
		insert_data(X, "END");

		std::cout << "Insert Y (write END to stop):" << std::endl;
		insert_data(Y, "END");

		std::cout << "Error on Y:" << std::endl;
		std::cin >> sigma_Y;

		real A = lst_sqrs_lin_intercept(X, Y);
		real B = lst_sqrs_lin_slope(X, Y);

		if(sigma_Y == 0) {
			sigma_Y = least_squares_linear_error(X, Y, A, B);
		}

		std::cout << "Covariance = " << sample_covariance(X, Y) << std::endl;
		std::cout << "Pearson's Correlation Coefficient = "
			<< sample_correlation_coefficient(X, Y) << std::endl;
		std::cout << "r-Squared = " << square(sample_correlation_coefficient(X, Y)) << std::endl;

		std::cout << "\nOrdinary Least Squares Linearization:" << std::endl;
		std::cout << "A = " << A << " +- " << least_squares_linear_sigma_A(X, Y, sigma_Y) << std::endl;
		std::cout << "B = " << B << " +- " << least_squares_linear_sigma_B(X, Y, sigma_Y) << std::endl;
		std::cout << "Linearization Error = " << least_squares_linear_error(X, Y, A, B) << std::endl;

	} else if(n == 3) {

		std::cout << "Insert X (write END to stop): " << std::endl;
		insert_data(X, "END");

		std::cout << "Insert Y (write END to stop): " << std::endl;
		insert_data(Y, "END");

		std::cout << "Insert Z (write END to stop): " << std::endl;
		insert_data(Z, "END");

		real A = lst_sqrs_weight_lin_intercept(X, Y, Z);
		real B = lst_sqrs_weight_lin_slope(X, Y, Z);

		std::cout << "Covariance = " << sample_covariance(X, Y) << std::endl;
		std::cout << "Pearson's Correlation Coefficient = "
			<< sample_correlation_coefficient(X, Y) << std::endl;
		std::cout << "r-Squared = " << square(sample_correlation_coefficient(X, Y)) << std::endl;

		std::cout << "\nOrdinary Least Squares Linearization:" << std::endl;
		std::cout << "A = " << A << " +- " << least_squares_linear_sigma_A(X, Y, sigma_Y) << std::endl;
		std::cout << "B = " << B << " +- " << least_squares_linear_sigma_B(X, Y, sigma_Y) << std::endl;
		std::cout << "Linearization Error = " << least_squares_linear_error(X, Y, A, B) << std::endl;

		std::cout << "Linearization Chi-Square = " << chi_square_linearization(X, Y, Z, A, B) << std::endl;
		
	} else {
		std::cout << "Input error" << std::endl;
		return 1;
	}

	return 0;
}
