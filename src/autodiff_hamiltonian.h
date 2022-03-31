
///
/// @file autodiff_hamiltonian.cpp Automatically simulate an Hamiltonian system from
/// its Hamiltonian function using automatical differentiation and
/// numerical integration
///

#include <fstream>
#include "uroboro.h"

using namespace umath;


dual harmonic_oscillator(vec<2, dual> eta) {

	const real omega = 1;
	const real m = 1;

	return square(eta[1]) / (2.0 * m) + square(eta[0]) * square(omega) * m / 2.0;
}


dual pendulum(vec<2, dual> eta) {

	const real l = 1;
	const real m = 1;
	const real g = 9.81;

	return square(eta[1]) / (2.0 * m) + (-cos(eta[0]) + 1) * m * g * l;
}


int main(int argc, char const *argv[]) {

	std::cout.precision(12);

	vec2 eta = {0, 1};

	vec2 eta_dt;
	const real dt = 0.0001;

	std::ofstream of("harmonic_oscillator.dat");

	for (int i = 0; i < 1000000; ++i) {

		eta_dt = mat2::symplectic() * gradient(harmonic_oscillator, eta);
		eta = eta + eta_dt * dt;
		
		of << eta[0] << " " << eta[1] << std::endl;
	}

	of.close();

	return 0;
}
