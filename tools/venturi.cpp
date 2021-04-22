#define UROBORO_INCLUDE_ALL
#include "./statistics.h"
#include "./utility.h"
#include <iostream>
#include <string>

using namespace uroboro;

real square(real x) {
    return x * x;
}

real portata(real m, real rho, real t) {
    return m / (rho * t);
}

real portata(real A, real v) {
    return A * v;
}

real ps_bernoulli(real pl, real rho, real As, real Al, real R) {
    return pl - 0.5 * rho * ((1 / square(As)) - (1/square(Al))) * square(R);
}

real ps_corr(real p1, real p2, real p3) {
    return p2 + (p1 - p3) / 2.0;
}


int main(int argc, char const *argv[]) {

	std::cout << "Insert each value and press Enter, write END to stop insertion" << std::endl;

	vec_buff X;
	std::string line;
	real value;

	while(true) {

		std::getline(std::cin, line);

		if(line == "END")
			break;

		try {
			value = std::stod(line);
		} catch(...) {
			std::cout << "Input conversion error" << std::endl;
		}

		X.emplace_back(value);

	}

	print_sample_stats(X);

	real p1 = X[0];
	real p2 = X[1];
	real p3 = X[2];
	real p4 = X[3];

	// Costanti
	real rho_acqua = 1000;
	real rho_aria = 1.225;
	real As = 0.0000452;
	real Al = 0.000199;

	// Variabili
	// real R = portata(m, rho_acqua, t);
	real R = 0.000016;

    std::cout << ps_bernoulli(p1, rho_aria, As, Al, R) << std::endl;
    std::cout << ps_corr(p1, p2, p3) << std::endl;
    std::cout << uroboro::abs(ps_bernoulli(p1, rho_aria, As, Al, R) - ps_corr(p1, p2, p3)) << std::endl;

	std::cout << "Press Enter to exit..." << std::endl;
	std::cin.get();

	return 0;
}
