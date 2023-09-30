
///
/// @file hamiltonian.cpp Automatically simulate an Hamiltonian system from
/// its Hamiltonian function using automatic differentiation and numerical integration.
/// This example may be compiled using 'make hamiltonian_simulation'
///

#include <fstream>
#include "theoretica.h"
using namespace th;


// Dimension of configuration space
constexpr unsigned int N = 3;

// Dimension of phase space
constexpr unsigned int M = 2 * N;

// M-dimensional symplectic matrix
const auto M_symplectic = mat<real, M, M>::symplectic();


// Hamiltonian for the N-dimensional harmonic oscillator (m = 1, omega = 1)
template<typename NumType>
NumType harmonic_oscillator(vec<M, NumType> eta) {

	return (eta * eta) / 2.0;
}


// Hamiltonian for the simple pendulum
template<typename NumType>
NumType pendulum(vec<2, NumType> eta) {

	return square(eta[1]) / 2.0 - cos(eta[0]);
}


// Differential vector field
vec<M> f(real t, vec<M> eta) {

	return M_symplectic * gradient<M>(harmonic_oscillator, eta);
}


int main(int argc, char const *argv[]) {

	// Initial conditions
	vec<M> s0 = {0, 1, 2, 1, 0.5, 0};

	// State of the system
	ode_state<M> state(s0);

	// Time step
	const real dt = 0.001;
	
	// Number of iterations
	const unsigned int iterations = 50000;

	// Output file
	std::ofstream file("hamiltonian.dat");

	for (unsigned int i = 0; i < iterations; ++i) {

		state = ode_rk4(f, state, dt);
		file << state << std::endl;
	}

	file.close();

	return 0;
}
