#ifndef THEORETICA_QUANTUM_H
#define THEORETICA_QUANTUM_H

#include "../complex/complex.h"
#include "../polynomial/ortho_polyn.h"


namespace theoretica {


	// Constants

	// Electron mass in kilograms
	constexpr real ELECTRON_MASS = 9.1093837015E-31;

	// Proton mass in kilograms
	constexpr real PROTON_MASS = 1.67262192369E-27;

	// Bohr radius in nanometers
	constexpr real BOHR_RADIUS = 5.29177210903;

	// Reduced Bohr radius in nanometers
	constexpr real BOHR_RADIUS_REDUC = BOHR_RADIUS * ((ELECTRON_MASS + PROTON_MASS) / (PROTON_MASS));


	/// One-dimensional complex wave function
	using wavefunc_1D = std::function<complex(real)>;
	
	/// Two-dimensional complex wave function
	using wavefunc_2D = std::function<complex(real, real)>;
	
	/// Three-dimensional complex wave function
	using wavefunc_3D = std::function<complex(real, real, real)>;


	/// Spherical harmonics normalized for quantum mechanics
	inline std::function<complex(real, real)> spherical_harmonic_qm(unsigned int l, int m) {

		if(abs(m) > l) {
			TH_MATH_ERROR("spherical_harmonic", m, IMPOSSIBLE_OPERATION);
			return [](real a, real b) { return complex(nan(), nan()); };
		}

		const auto Leg = assoc_legendre_polynomial(l, m);
		const real K = sqrt(
			(2 * ((real) l) + 1) * (fact(l - m) / (real) fact(l + m)) / (4 * PI));

		if(m != 0) {

			return [=](real theta, real phi) {
				return K * Leg(cos(theta)) * exp(complex::i() * m * phi);
			};

		} else {

			return [=](real theta, real phi) {
				return K * Leg(cos(theta));
			};
		}
	}


	/// Complex wave function of the hydrogen atom
	/// @param n The principal quantum number
	/// @param l The azimuthal  quantum number
	/// @param m The magnetic quantum number
	/// @return A lambda function representing the wave function of the hydrogen atom
	inline wavefunc_3D wavefunc_hydrogen(unsigned int n, unsigned int l, int m) {

		// Check quantum numbers
		if(n == 0) {
			TH_MATH_ERROR("wavefunc_hydrogen", n, IMPOSSIBLE_OPERATION);
			return [](real a, real b, real c) { return complex(nan(), nan()); };
		}

		if(l >= n) {
			TH_MATH_ERROR("wavefunc_hydrogen", l, IMPOSSIBLE_OPERATION);
			return [](real a, real b, real c) { return complex(nan(), nan()); };
		}

		if(abs(m) > l) {
			TH_MATH_ERROR("wavefunc_hydrogen", m, IMPOSSIBLE_OPERATION);
			return [](real a, real b, real c) { return complex(nan(), nan()); };
		}


		// Normalization constant
		const real N = sqrt(pow(2.0 / (n * BOHR_RADIUS_REDUC), 3)
			* fact(n - l - 1) / (2 * n * fact(n + l)));

		const polynomial<real> Lag = general_laguerre_polynomial(2 * l + 1, n - l - 1);
		const auto Y = spherical_harmonic_qm(l, m);

		return [=](real r, real theta, real phi) {

			const real rho = (2 * r) / (n * BOHR_RADIUS_REDUC);
			return N * pow(rho, l) * Lag(rho) * exp(-rho / 2) * Y(theta, phi);
		};
	}


}


#endif
