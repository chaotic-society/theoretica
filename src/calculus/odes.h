
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

			static_assert(N > 0, "N cannot be zero.");

			real t;
			vec<N> y;

			ode_state() : t(0), y(0) {}

			ode_state(real t, vec<N> y) : t(t), y(y) {}

			ode_state(vec<N> y) : t(0), y(y) {}

			inline ode_state<N>& operator=(const std::array<real, 4>& v) {
				t = v[0];
				y[0] = v[1];
				y[1] = v[2];
				y[2] = v[3];
				return *this;
			}

#ifndef THEORETICA_NO_PRINT

			/// Convert the ODE state to string representation
			inline std::string to_string(const std::string& separator = " ") const {

				std::stringstream res;
				res << t << separator;

				for (unsigned int i = 0; i < N; ++i) {

					res << y.data[i];

					if(i != N - 1)
						res << separator;
				}

				return res.str();
			}


			/// Stream the ODE state in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const ode_state<N>& obj) {
				return out << obj.to_string();
			}
#endif
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

			ode_state(real y) : t(0), y(y) {}

			inline ode_state<1>& operator=(const std::array<real, 2>& v) {
				t = v[0];
				y = v[1];
				return *this;
			}

#ifndef THEORETICA_NO_PRINT

			/// Convert the ODE state to string representation
			inline std::string to_string(const std::string& separator = " ") const {

				std::stringstream res;
				res << t << separator << y;
				return res.str();
			}


			/// Stream the ODE state in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const ode_state<1>& obj) {
				return out << obj.to_string();
			}
#endif
	};


	/// Integrate numerically a differential equation in one unknown
	/// using Euler's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	inline ode_state<1> ode_euler(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(s.t + h, s.y + h * f(s.t, s.y));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Euler's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_euler(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		return ode_state<N>(s.t + h, s.y + h * f(s.t, s.y));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using the midpoint method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	inline ode_state<1> ode_midpoint(real(*f)(real, real), ode_state<1> s, real h) {

		return ode_state<1>(
			s.t + h,
			s.y + h * f(s.t + h / 2.0, s.y + f(s.t, s.y) * h / 2.0));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using the midpoint method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_midpoint(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		return ode_state<N>(
			s.t + h,
			s.y + h * f(s.t + h / 2.0, s.y + f(s.t, s.y) * h / 2.0));
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Heun's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	inline ode_state<1> ode_heun(real(*f)(real, real), ode_state<1> s, real h) {

		const real y_p = s.y + h * f(s.t, s.y);
		const real t_new = s.t + h;

		return ode_state<1>(t_new, s.y + (f(s.t, s.y) + f(t_new, y_p)) * h / 2.0);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Heun's method.
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_heun(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> y_p = s.y + h * f(s.t, s.y);
		const real t_new = s.t + h;

		return ode_state<N>(t_new, s.y + (f(s.t, s.y) + f(t_new, y_p)) * h / 2.0);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Runge-Kutta's method of second order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	inline ode_state<1> ode_rk2(real(*f)(real, real), ode_state<1> s, real h) {

		const real k1 = f(s.t, s.y);
		const real k2 = f(s.t + (h / 2.0), s.y + k1 * (h / 2.0));

		return ode_state<1>(s.t + h, s.y + k2 * h);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using Runge-Kutta's method of second order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_rk2(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> k1 = f(s.t, s.y);
		const vec<N> k2 = f(s.t + (h / 2.0), s.y + k1 * (h / 2.0));

		return ode_state<N>(s.t + h, s.y + k2 * h);
	}


	/// Integrate numerically a differential equation in one unknown
	/// using Runge-Kutta's method of fourth order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	inline ode_state<1> ode_rk4(real(*f)(real, real), ode_state<1> s, real h) {

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
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_rk4(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

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
	/// @param h Step size
	inline ode_state<1> ode_k38(real(*f)(real, real), ode_state<1> s, real h) {

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
	/// @param h Step size
	template<unsigned int N>
	inline ode_state<N> ode_k38(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

		const vec<N> k1 = f(s.t, s.y);
		const vec<N> k2 = f(s.t + (h / 3.0), s.y + k1 * (h / 3.0));
		const vec<N> k3 = f(s.t + (h * 2 / 3.0), s.y + h * (-k1 / 3.0 + k2));
		const vec<N> k4 = f(s.t + h, s.y + h * (k1 - k2 + k1));

		return ode_state<N>(s.t + h, s.y + (k1 + 3 * k2 + 3 * k3 + k4) * h / 8.0);
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using the Dormand-Prince method of sixth order.
	///
	/// The Dormand-Prince method is applied with the initial step size and
	/// further iterations are computed until a small enough step size
	/// with an estimated error smaller than the tolerance is found,
	/// halving the step size at each iteration. If the estimated error
	/// is less than half the tolerance, the step size is doubled to speed
	/// up calculations. 
	///
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	/// @param tolerance The tolerance between fourth and fifth order solutions
	// inline ode_state<1> integrate_ode_rkdp(real_function f, ode_state<N> s,
	// 	real& h, real tolerance = THEORETICA_RKDP_TOL) {

	// 	const real k1 = f(s.t, s.y);
	// 	const real k2 = f(s.t + h / 5.0,
	// 		s.y + h * (k1 * 1.0 / 5.0));

	// 	const real k3 = f(s.t + h * 3.0 / 10.0,
	// 		s.y + h * (k1 * 3.0 / 40.0 + k2 * 9.0 / 40.0));

	// 	const real k4 = f(s.t + h * 4.0 / 5.0,
	// 		s.y + h * (k1 * 44.0 / 45.0 - k2 * 56.0 / 15.0 + k3 * 32.0 / 9.0));

	// 	const real k5 = f(s.t + h * 8.0 / 9.0,
	// 		s.y + h * (k1 * 19372.0/6561.0 - k2 * 25360.0/2187.0 + k3 * 64448.0/6561.0 - k4 * 212.0/729.0));

	// 	const real k6 = f(s.t + h,
	// 		s.y + h * (k1 * 9017.0/3168.0 - k2 * 355.0/33.0 + k3 * 46732.0/5247.0 + k4 * 49.0/176.0 - k5 * 5103.0/18656.0));

	// 	const real k7 = f(s.t + h,
	// 		s.y + h * (k1 * 35.0/384 + k3 * 500.0/1113.0 + k4 * 125.0/192.0 - k5 * 2187.0/6784 + k6 * 11.0/84.0));

	// 	real y_4 = k7;
	// 	real y_best = s.y + h * (
	// 		  k1 * 5179.0/57600.0	+ k3 * 7571.0/16695.0
	// 		+ k4 * 393.0/640.0		- k5 * 92097.0 / 339200.0
	// 		+ k6 * 187.0/2100.0		+ k7 / 40.0);

	// 	return ode_state<1>(s.t + h, s.y + h * (
	// 		  k1 * 5179.0/57600.0 + k3 * 7571.0/16695.0
	// 		+ k4 * 393.0/640.0 - k5 * 92097.0 / 339200.0
	// 		+ k6 * 187.0/2100.0 + k7 / 40.0));
	// }


	/// Integrate numerically a differential equation in N unknowns
	/// using the Dormand-Prince method of sixth order
	/// @param f A function which computes the derivative of the state
	/// @param s The current state of the ODE
	/// @param h Step size
	// template<unsigned int N>
	// inline ode_state<N> ode_rkdp(vec<N>(*f)(real, vec<N>), ode_state<N> s, real h) {

	// 	const vec<N> k1 = f(s.t, s.y);
	// 	const vec<N> k2 = f(s.t + h / 5.0,
	// 		s.y + h * (k1 * 1.0 / 5.0));

	// 	const vec<N> k3 = f(s.t + h * 3.0 / 10.0,
	// 		s.y + h * (k1 * 3.0 / 40.0 + k2 * 9.0 / 40.0));

	// 	const vec<N> k4 = f(s.t + h * 4.0 / 5.0,
	// 		s.y + h * (k1 * 44.0 / 45.0 - k2 * 56.0 / 15.0 + k3 * 32.0 / 9.0));

	// 	const vec<N> k5 = f(s.t + h * 8.0 / 9.0,
	// 		s.y + h * (k1 * 19372.0/6561.0 - k2 * 25360.0/2187.0 + k3 * 64448.0/6561.0 - k4 * 212.0/729.0));

	// 	const vec<N> k6 = f(s.t + h,
	// 		s.y + h * (k1 * 9017.0/3168.0 - k2 * 355.0/33.0 + k3 * 46732.0/5247.0 + k4 * 49.0/176.0 - k5 * 5103.0/18656.0));

	// 	const vec<N> k7 = f(s.t + h,
	// 		s.y + h * (k1 * 35.0/384 + k3 * 500.0/1113.0 + k4 * 125.0/192.0 - k5 * 2187.0/6784 + k6 * 11.0/84.0));

	// 	return ode_state<N>(s.t + h, s.y + h * (
	// 		  k1 * 5179.0/57600.0 + k3 * 7571.0/16695.0
	// 		+ k4 * 393.0/640.0 - k5 * 92097.0 / 339200.0
	// 		+ k6 * 187.0/2100.0 + k7 / 40.0));
	// }


	/// Integrate numerically a differential equation in 1 unknown
	/// using the Adams-Bashforth linear multistep method of the second order
	inline ode_state<1> ode_adams(real(*f)(real, real),
		ode_state<1> s0, ode_state<1> s1, real h) {

		return ode_state<1>(
			s1.t + h,
			s1.y + h * (3. * f(s1.t, s1.y) / 2. - f(s0.t, s0.y) / 2.));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using the Adams-Bashforth linear multistep method of the second order
	template<unsigned int N>
	inline ode_state<N> ode_adams(vec<N>(*f)(real, vec<N>),
		ode_state<N> s0, ode_state<N> s1, real h) {

		return ode_state<N>(
			s1.t + h,
			s1.y + h * (3. * f(s1.t, s1.y) / 2. - f(s0.t, s0.y) / 2.));
	}


	/// Integrate numerically a differential equation in 1 unknown
	/// using the Adams-Bashforth linear multistep method of third order
	inline ode_state<1> ode_adams3(real(*f)(real, real),
		ode_state<1> s0, ode_state<1> s1, ode_state<1> s2, real h) {

		return ode_state<1>(
			s2.t + h,
			s2.y + h * (23. / 12. * f(s2.t, s2.y)
				- 4. / 3.0 * f(s1.t, s1.y)
				+ 3. / 8. * f(s0.t, s0.y)));
	}


	/// Integrate numerically a differential equation in N unknowns
	/// using the Adams-Bashforth linear multistep method of third order
	template<unsigned int N>
	inline ode_state<N> ode_adams3(vec<N>(*f)(real, vec<N>),
		ode_state<N> s0, ode_state<N> s1, ode_state<N> s2, real h) {

		return ode_state<N>(
			s2.t + h,
			s2.y + h *	(23. / 12. * f(s2.t, s2.y)
						- 4. / 3.0 * f(s1.t, s1.y)
						+ 3. / 8. * f(s0.t, s0.y)));
	}

}

#endif
