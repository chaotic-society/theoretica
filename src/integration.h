#ifndef UROBORO_INTEGRATION_H
#define UROBORO_INTEGRATION_H

#include "./constants.h"
#include "./function.h"
#include "./polynomial.h"


namespace uroboro {


	// Integrate a polynomial
	polynomial integrate_polynomial(polynomial p) {

		polynomial Dp;
		Dp.coeff.push_back(0);

		for (int i = 0; i < p.size(); ++i) {
			Dp.coeff.push_back(p.get(i) / (real) (i + 1));
		}

		return Dp;
	}


	// Approximate the definite integral of an arbitrary function
	// using the midpoint method
	real approx_integral_midpoint(real_function f, real a, real b, unsigned int steps = 100) {
		
		real dx = (b - a) / steps;
		real res = 0;

		for (int i = 0; i < steps; ++i) {
			res += f(a + (i + 0.5) * dx);
		}

		return res * dx;
	}


	// Approximate the definite integral of an arbitrary function
	// using the trapezoid method
	real approx_integral_trapezoid(real_function f, real a, real b, unsigned int steps = 100) {
		
		real dx = (b - a) / steps;
		real res = 0;

		res += 0.5 * f(a);

		for (int i = 1; i < steps; ++i) {
			res += f(a + i * dx);
		}

		res += 0.5 * f(b);

		return res * dx;
	}


	// Approximate the definite integral of an arbitrary function
	// using Simpson's method
	real approx_integral_simpson(real_function f, real a, real b, unsigned int steps = 100) {
		
		real dx = (b - a) / (real) steps;
		real res = 0;

		res += f(a) + f(b);

		for (int i = 1; i < steps; ++i) {

			if(i % 2 == 0)
				res += 2.0 * f(a + i * dx);
			else
				res += 4.0 * f(a + i * dx);
			
		}

		return res * dx / 3.0;
	}


	// TO-DO Romberg integration

	// TO-DO Monte Carlo integration


	// Runge-Kutta integration of 4th order
	namespace RK4 {

		// A kinematic state (position, velocity)
		struct kinematic_state {

			real x;
			real v;

			kinematic_state() : x(0), v(0) {}
			kinematic_state(real x, real v) : x(x), v(v) {}
		};


		// A kinematic state derivative (dx, dv)
		struct kinematic_deriv {

			real dx;
			real dv;

			kinematic_deriv() : dx(0), dv(0) {}
			kinematic_deriv(real dx, real dv) : dx(dx), dv(dv) {}
		};


		// Function type for acceleration functions
		using accel_function = real(*)(const kinematic_state&, real);


		// Eval a single kinematic state with its derivative
		inline kinematic_deriv eval(
			const kinematic_state& prec,
			real t, real dt,
			const kinematic_deriv& deriv,
			accel_function accel) {

			kinematic_state state = kinematic_state(
				prec.x + deriv.dx * dt,
				prec.v + deriv.dv * dt);

			kinematic_deriv res = kinematic_deriv(state.v, accel(state, t + dt));

			return res;
		}


		// Runge-Kutta integration of 4th order
		inline void integrate(kinematic_state s, real t, real dt, accel_function accel) {

			kinematic_deriv a, b, c, d;

			a = eval(s, t, 0, kinematic_deriv(), accel);
			b = eval(s, t, dt * 0.5, a, accel);
			c = eval(s, t, dt * 0.5, b, accel);
			d = eval(s, t, dt, c, accel);

			// Approximate integration using 4 different points
			real dxdt = 1.0f / 6.0f * (a.dx + 2.0f * (b.dx + c.dx) + d.dx);
			real dvdt = 1.0f / 6.0f * (a.dv + 2.0f * (b.dv + c.dv) + d.dv);

			// Update state
			s.x = s.x + dxdt * dt;
			s.v = s.v + dvdt * dt;
		}

	}

}


#endif
