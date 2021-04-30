#ifndef UROBORO_INTEGRATION_H
#define UROBORO_INTEGRATION_H

#include "./constants.h"
#include "./function.h"



namespace uroboro {


	real integrate_approx(real_function f, real a, real b, unsigned int prec = 1000) {

		// TO-DO Simpson Cavalieri approximation

		return 0;
	}


	// Runge-Kutta integration of 4th order

	struct kinematic_state {

		real x;
		real v;

		kinematic_state() : x(0), v(0) {}
		kinematic_state(real x, real v) : x(x), v(v) {}

	};


	struct kinematic_deriv {

		real dx;
		real dv;

		kinematic_deriv() : dx(0), dv(0) {}
		kinematic_deriv(real dx, real dv) : dx(dx), dv(dv) {}

	};


	using acceleration_function = real(*)(const kinematic_state&, real);


	inline kinematic_deriv eval(
		const kinematic_state& prec,
		real t, real dt,
		const kinematic_deriv& deriv,
		acceleration_function accel) {

		kinematic_state state = kinematic_state(
			prec.x + deriv.dx * dt,
			prec.v + deriv.dv * dt);

		kinematic_deriv res = kinematic_deriv(state.v, accel(state, t + dt));

		return res;
	}


	inline void integrate_rk4(kinematic_state s, real t, real dt, acceleration_function accel) {

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


#endif
