
///
/// @file odes.h ODE integration methods
///

#ifndef THEORETICA_ODES_H
#define THEORETICA_ODES_H

#include "../core/constants.h"
#include <tuple>


namespace theoretica {


	/// @class ode_state
	/// The current state of an ODE integration
	/// for an N dimensional differential equation.
	template<unsigned int N>
	class ode_state {
		public:
			real t;
			vec<N> y;

			ode_state(real t, vec<N> y) {
				this->t = t;
				this->y = y;
			}
	};


	/// @class ode_state
	/// The current state of an ODE integration
	/// for a differential equation in one unknown.
	template<>
	class ode_state<1> {
	public:
			real t;
			real y;

			ode_state(real t, real y) {
				this->t = t;
				this->y = y;
			}
	};


	/// Integrate numerically a differential equation in one unknown
	/// using Euler's method.
	ode_state<1> integrate_ode_euler(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(s.t + h, s.y + h * f(s.t, s.y));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using the midpoint method.
	ode_state<1> integrate_ode_midpoint(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(
			s.t + h,
			s.y + h * f(s.t + h / 2.0, s.y + f(s.t, s.y) * h / 2.0));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Heun's method.
	ode_state<1> integrate_ode_heun(real(*f)(real, real), ode_state<1> s, real h) {

		const real y_p = s.y + h * f(s.y, s.y);
		const real t_new = s.t + h;

		return ode_state<1>(t_new, s.y + (f(s.t, s.y) + f(t_new, y_p)) * h / 2.0);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Runge-Kutta's method of fourth order
	ode_state<1> integrate_ode_rk4(real(*f)(real, real), ode_state<1> s, real h) {

		const real k1 = f(s.t, s.y);
		const real k2 = f(s.t + (h / 2.0), s.y + k1 * (h / 2.0));
		const real k3 = f(s.t + (h / 2.0), s.y + k2 * (h / 2.0));
		const real k4 = f(s.t + h, s.y + k3 * h / 2.0);

		return ode_state<1>(s.t + h, s.y + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6.0);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Kutta's 3/8 rule method
	ode_state<1> integrate_ode_k38(real(*f)(real, real), ode_state<1> s, real h) {

		const real k1 = f(s.t, s.y);
		const real k2 = f(s.t + (h / 3.0), s.y + k1 * (h / 3.0));
		const real k3 = f(s.t + (h * 2 / 3.0), s.y + h * (-k1 / 3.0 + k2));
		const real k4 = f(s.t + h, s.y + h * (k1 - k2 + k1));

		return ode_state<1>(s.t + h, s.y + (k1 + 3 * k2 + 3 * k3 + k4) * h / 8.0);
	}



}

#endif
