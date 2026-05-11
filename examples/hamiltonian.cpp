
///
/// @file hamiltonian.cpp Simulate an Hamiltonian system from
/// its Hamiltonian function using automatic differentiation and numerical integration.
/// This example may be compiled using 'make example_hamiltonian'
///

#include <fstream>
#include "theoretica.h"
using namespace th;
using namespace autodiff;


/// Given the Hamiltonian of a conservative system, computes its Hamiltonian flow of the system.
template<unsigned int N>
inline auto hamiltonian_flow(std::function<dreal_t<N>(const dvec_t<N>&)> H) {

	using Vector = vec<real, N>;

	return [H](real t, const Vector& z) {

		Vector flow;
		flow.resize(z.size());

		// Compute gradient of the Hamiltonian
		auto grad = autodiff::gradient(H, z);

		// Symplectic structure
		const size_t half_size = flow.size() / 2;
		size_t i = 0;

		for (; i < half_size; i++) {
			flow[i] = grad[i + half_size];
			flow[i + half_size] = -grad[i];
		}

		return flow;
	};
}


// Hamiltonian for the 2D harmonic oscillator (m = 1, omega = 1)
dreal4 harmonic_oscillator(dvec4 eta) {

	return (eta * eta) / 2.0;
}


// Hamiltonian for the simple pendulum
dreal2 pendulum(dvec2 eta) {

	return square(eta[1]) / 2.0 - cos(eta[0]);
}


int main() {

	// Initial conditions
	vec4 s0 = {
		0.0, 1.0,
		1.0, 0.5,
	};

	// Time step
	const real dt = 0.001;

	// Final time
	const real tf = 50.0;

	// Output file
	std::ofstream file ("hamiltonian.dat");

	// Compute the trajectory using Runge-Kutta's method.
	auto solution = ode::solve_rk4(
		hamiltonian_flow<4>(harmonic_oscillator),
		s0, 0.0, tf, dt
	);

	file << solution;
}
