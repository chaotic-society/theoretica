
///
/// @file hamiltonian.cpp Automatically simulate an Hamiltonian system from
/// its Hamiltonian function using automatic differentiation and numerical integration.
/// This example may be compiled using 'make hamiltonian_simulation'
///

#include <fstream>
#include "theoretica.h"
using namespace th;
using namespace autodiff;


// Dimension of configuration space
constexpr unsigned int N = 3;

// Dimension of phase space
constexpr unsigned int M = 2 * N;

// MxM symplectic matrix
const auto M_symplectic = algebra::symplectic<mat<real, M, M>>();


// Hamiltonian for the N-dimensional harmonic oscillator (m = 1, omega = 1)
dreal_t<M> harmonic_oscillator(dvec_t<M> eta) {

	return (eta * eta) / 2.0;
}


// Hamiltonian for the simple pendulum
template<typename NumType>
dreal2 pendulum(dvec2 eta) {

	return square(eta[1]) / 2.0 - cos(eta[0]);
}


// Differential vector field
vec<real, M> f(real t, vec<real, M> eta) {

	return M_symplectic * autodiff::gradient(harmonic_oscillator, eta);
}


int main() {

	// Initial conditions
	vec<real, M> s0 = {
		0.0, 1.0, 2.0,
		1.0, 0.5, 0.0
	};

	// Time step
	const real dt = 0.001;

	// Final time
	const real tf = 50.0;

	// Output file
	std::ofstream file("hamiltonian.dat");

	// Compute the trajectory using Runge-Kutta's method.
	auto solution = ode::solve_rk4(f, s0, 0.0, tf, dt);
}
