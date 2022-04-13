
///
/// @file error_propagation.h Automatic error propagation
///

#include "uroboro.h"
using namespace umath;


// Example function to compute error on
template<typename NumType>
NumType theta(vec<2, NumType> v) {

	const NumType r = v[0];
	const NumType d = v[1];

	return umath::atan(r / (d - 1));
}


int main(int argc, char const *argv[]) {

	std::cout.precision(12);

	real r, d;
	real delta_r, delta_d_lf;

	std::cout << "r = ";
	std::cin >> r;

	std::cout << "delta_r = ";
	std::cin >> delta_r;

	std::cout << "d = ";
	std::cin >> d;

	std::cout << "delta_d_lf = ";
	std::cin >> delta_d_lf;

	vec2 data = {r, d};
	vec2 err = {delta_r, delta_d_lf};

	auto grad = gradient(theta, data);
	real err_tot = 0;

	// Sum under quadrature
	for (int i = 0; i < data.size(); ++i)
		err_tot += square(grad.at(i)) * square(err.at(i));

	err_tot = umath::sqrt(err_tot);

	// Linear sum
	// for (int i = 0; i < data.size(); ++i)
	// 	err_tot += grad.at(i) * err.at(i);

	std::cout << "Total error:\n" << err_tot << std::endl;

	return 0;
}
