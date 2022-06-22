
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
	struct ode_state {

			real t;
			vec<N> y;

			ode_state() : t(0), y(0) {}

			ode_state(real t, vec<N> y) : t(t), y(y) {}
	};


	/// @class ode_state
	/// The current state of an ODE integration
	/// for a differential equation in one unknown.
	template<>
	struct ode_state<1> {

			real t;
			real y;

			ode_state() : t(0), y(0) {}

			ode_state(real t, real y) : t(t), y(y) {}
	};


	/// Integrate numerically a differential equation in one unknown
	/// using Euler's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	ode_state<1> integrate_ode_euler(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(s.t + h, s.y + h * f(s.t, s.y));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Euler's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	template<unsigned int N>
	ode_state<N> integrate_ode_euler(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		return ode_state<N>(s.t + h, s.y + h * f(s.t, s.y));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using the midpoint method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	ode_state<1> integrate_ode_midpoint(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(
			s.t + h,
			s.y + h * f(s.t + h / 2.0, s.y + f(s.t, s.y) * h / 2.0));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using the midpoint method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	template<unsigned int N>
	ode_state<N> integrate_ode_midpoint(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		return ode_state<N>(
			s.t + h,
			s.y + h * f(s.t + h / 2.0, s.y + f(s.t, s.y) * h / 2.0));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Heun's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	ode_state<1> integrate_ode_heun(real(*f)(real, real), ode_state<1> s, real h) {

		const real y_p = s.y + h * f(s.t, s.y);
		const real t_new = s.t + h;

		return ode_state<1>(t_new, s.y + (f(s.t, s.y) + f(t_new, y_p)) * h / 2.0);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Heun's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	template<unsigned int N>
	ode_state<N> integrate_ode_heun(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> y_p = s.y + h * f(s.t, s.y);
		const real t_new = s.t + h;

		return ode_state<N>(t_new, s.y + (f(s.t, s.y) + f(t_new, y_p)) * h / 2.0);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Runge-Kutta's method of fourth order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	ode_state<1> integrate_ode_rk4(real(*f)(real, real), ode_state<1> s, real h) {

		const real k1 = f(s.t, s.y);
		const real k2 = f(s.t + (h / 2.0), s.y + k1 * (h / 2.0));
		const real k3 = f(s.t + (h / 2.0), s.y + k2 * (h / 2.0));
		const real k4 = f(s.t + h, s.y + k3 * h / 2.0);

		return ode_state<1>(s.t + h, s.y + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6.0);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Runge-Kutta's method of fourth order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	template<unsigned int N>
	ode_state<N> integrate_ode_rk4(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> k1 = f(s.t, s.y);
		const vec<N> k2 = f(s.t + (h / 2.0), s.y + k1 * (h / 2.0));
		const vec<N> k3 = f(s.t + (h / 2.0), s.y + k2 * (h / 2.0));
		const vec<N> k4 = f(s.t + h, s.y + k3 * h / 2.0);

		return ode_state<N>(s.t + h, s.y + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6.0);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Kutta's 3/8 rule method
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	ode_state<1> integrate_ode_k38(real(*f)(real, real), ode_state<1> s, real h) {

		const real k1 = f(s.t, s.y);
		const real k2 = f(s.t + (h / 3.0), s.y + k1 * (h / 3.0));
		const real k3 = f(s.t + (h * 2 / 3.0), s.y + h * (-k1 / 3.0 + k2));
		const real k4 = f(s.t + h, s.y + h * (k1 - k2 + k1));

		return ode_state<1>(s.t + h, s.y + (k1 + 3 * k2 + 3 * k3 + k4) * h / 8.0);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Kutta's 3/8 rule method
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h The real step size
	template<unsigned int N>
	ode_state<N> integrate_ode_k38(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> k1 = f(s.t, s.y);
		const vec<N> k2 = f(s.t + (h / 3.0), s.y + k1 * (h / 3.0));
		const vec<N> k3 = f(s.t + (h * 2 / 3.0), s.y + h * (-k1 / 3.0 + k2));
		const vec<N> k4 = f(s.t + h, s.y + h * (k1 - k2 + k1));

		return ode_state<N>(s.t + h, s.y + (k1 + 3 * k2 + 3 * k3 + k4) * h / 8.0);
	}


	/// Integrate numerically a differential equation in 1 unknown
	// using the Adams-Bashforth linear multistep method of the second order
	ode_state<1> integrate_ode_adams(real(*f)(real, real),
		ode_state<1> s0, ode_state<1> s1, real h) {

		return ode_state<1>(
			s1.t + h,
			s1.y + h * (3. * f(s1.t, s1.y) / 2. - f(s0.t, s0.y) / 2.));
	}


	/// Integrate numerically a differential equation in N unknowns
	// using the Adams-Bashforth linear multistep method of the second order
	template<unsigned int N>
	ode_state<N> integrate_ode_adams(vec<N>(*f)(real, vec<N>),
		ode_state<N> s0, ode_state<N> s1, real h) {

		return ode_state<N>(
			s1.t + h,
			s1.y + h * (3. * f(s1.t, s1.y) / 2. - f(s0.t, s0.y) / 2.));
	}


	/// Integrate numerically a differential equation in 1 unknown
	// using the Adams-Bashforth linear multistep method of third order
	ode_state<1> integrate_ode_adams3(real(*f)(real, real),
		ode_state<1> s0, ode_state<1> s1, ode_state<1> s2, real h) {

		return ode_state<1>(
			s2.t + h,
			s2.y + h * (23. / 12. * f(s2.t, s2.y)
				- 4. / 3.0 * f(s1.t, s1.y)
				+ 3. / 8. * f(s0.t, s0.y)));
	}


	/// Integrate numerically a differential equation in N unknowns
	// using the Adams-Bashforth linear multistep method of third order
	template<unsigned int N>
	ode_state<N> integrate_ode_adams3(vec<N>(*f)(real, vec<N>),
		ode_state<N> s0, ode_state<N> s1, ode_state<N> s2, real h) {

		return ode_state<N>(
			s2.t + h,
			s2.y + h *	(23. / 12. * f(s2.t, s2.y)
						- 4. / 3.0 * f(s1.t, s1.y)
						+ 3. / 8. * f(s0.t, s0.y)));
	}

}

#endif
