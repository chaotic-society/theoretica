
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

// MxM symplectic matrix
const auto M_symplectic = mat<real, M, M>::symplectic();


// Hamiltonian for the N-dimensional harmonic oscillator (m = 1, omega = 1)
template<typename NumType>
NumType harmonic_oscillator(vec<NumType, M> eta) {

	return (eta * eta) / 2.0;
}


// Hamiltonian for the simple pendulum
template<typename NumType>
NumType pendulum(vec<NumType, 2> eta) {

	return square(eta[1]) / 2.0 - cos(eta[0]);
}


// Differential vector field
vec<real, M> f(real t, vec<real, M> eta) {

	return M_symplectic * autodiff::gradient<M>(harmonic_oscillator, eta);
}


int main() {

	// Initial conditions
	vec<real, M> s0 = {0, 1, 2, 1, 0.5, 0};

	// State of the system
	ode_state<M> state(s0);

	// Time step
	const real dt = 0.001;
	
	// Number of iterations
	const unsigned int iterations = 50000;

	// Output file
	std::ofstream file("hamiltonian.dat");

	// Compute the trajectory using Runge-Kutta's method.
	for (unsigned int i = 0; i < iterations; ++i) {

		state = ode_rk4(f, state, dt);
		file << state << std::endl;
	}
	
}
