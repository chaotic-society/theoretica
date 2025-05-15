
///
/// @file ode.h Numerical methods for ordinary differential equations.
///

#ifndef THEORETICA_ODE_H
#define THEORETICA_ODE_H

#include "../algebra/vec.h"


namespace theoretica {

	/// @namespace theoretica::ode Numerical methods for ordinary differential equations.
	namespace ode {


		/// @class ode_solution_t
		/// Data structure holding the numerical solution of a discretized ODE,
		/// where the vector \f$\vec t\f$ represents the discrete time points
		/// (independent variable) and the vector \f$\vec x\f$ the discrete solution.
		template<typename Vector = vec<real>>
		struct ode_solution_t {

			/// A vector of the time values (independent variable).
			vec<real> t;

			/// A vector of the phase space values (solution).
			vec<Vector> x;


			/// Default constructor
			ode_solution_t() {}


			/// Prepare the structure for integration by
			/// specifying the number of total steps and the initial
			/// conditions and time.
			ode_solution_t(size_t steps, const Vector& x0, real t0) {
					
				t.resize(steps);
				x.resize(steps);

				t[0] = t0;
				x[0] = x0;
			}


#ifndef THEORETICA_NO_PRINT			

			/// Convert the ODE solution to string representation
			inline std::string to_string(const std::string& separator = " ") const {

				if (t.size() != x.size()) {
					TH_MATH_ERROR("ode_solution_t::to_string", t.size(), INVALID_ARGUMENT);
					return "";
				}

				std::stringstream res;

				for (unsigned int i = 0; i < t.size(); ++i) {
					
					res << t[i] << separator;

					for (unsigned int j = 0; j < x[i].size(); ++j) {

						res << x[i][j];

						if(j != x[i].size() - 1)
							res << separator;
						else
							res << "\n";
					}
				}

				return res.str();
			}


			/// Convert the ODE solution to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the ODE solution in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(
				std::ostream& out, const ode_solution_t<Vector>& obj) {
				return out << obj.to_string();
			}
#endif

		};


		/// The solution of an ODE with any number of variables.
		using ode_solution = ode_solution_t<vec<real>>;


		/// The solution of an ODE in 1 variable.
		using ode_solution1d = ode_solution_t<real>;


		/// The solution of an ODE in 2 variables.
		using ode_solution2d = ode_solution_t<vec2>;


		/// The solution of an ODE in 3 variables.
		using ode_solution3d = ode_solution_t<vec3>;


		/// The solution of an ODE in 4 variables.
		using ode_solution4d = ode_solution_t<vec4>;


		/// A function representing a system of differential equations,
		/// taking as input the time (independent variable) and the current
		/// value of the variables (dependent variables), returning
		/// the time derivatives of each variable, such as \f$f\f$ in 
		/// \f$\dot \vec x = f(t, \vec x)\f$.
		template<typename Vector>
		using ode_function = std::function<Vector(real, const Vector&)>;


		// Steppers
		// (functions which compute one iteration of a method)


		/// Compute one step of Euler's method for ordinary differential equations.
		/// This function is used in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_euler(OdeFunction f, const Vector& x, real t, real h = 0.0001) {

			return x + h * f(t, x);
		}


		/// Compute one step of the midpoint method for ordinary differential equations.
		/// This function is used in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_midpoint(OdeFunction f, const Vector& x, real t, real h = 0.0001) {
			
			return x + h * f(t + h / 2.0, x + f(t, x) * h / 2.0);
		}


		/// Compute one step of Heun's method for ordinary differential equations.
		/// This function is used in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_heun(OdeFunction f, const Vector& x, real t, real h = 0.001) {
			
			const Vector x_p = x + h * f(t, x);
			const real t_new = t + h;

			return x + (f(t, x) + f(t_new, x_p)) * (h / 2.0);
		}


		/// Compute one step of the Runge-Kutta method of 2nd order for
		/// ordinary differential equations. This function is used in solvers
		/// to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_rk2(OdeFunction f, const Vector& x, real t, real h = 0.001) {
			
			const Vector k1 = f(t, x);
			const Vector k2 = f(t + (h / 2.0), x + k1 * (h / 2.0));

			return x + k2 * h;
		}


		/// Compute one step of the Runge-Kutta method of 4th order for
		/// ordinary differential equations. This function is used in solvers
		/// to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_rk4(OdeFunction f, const Vector& x, real t, real h = 0.01) {

			const real half = h / 2.0;
			
			const Vector k1 = f(t, x);
			const Vector k2 = f(t + half, x + k1 * half);
			const Vector k3 = f(t + half, x + k2 * half);
			const Vector k4 = f(t + h, x + k3 * h);

			return x + (k1 + k2 * 2.0 + k3 * 2.0 + k4) * (h / 6.0);
		}


		/// Compute one step of Kutta's 3/8 rule method for
		/// ordinary differential equations. This function is used
		/// in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_k38(OdeFunction f, const Vector& x, real t, real h = 0.001) {
			
			const Vector k1 = f(t, x);
			const Vector k2 = f(t + (h / 3.0), x + k1 * (h / 3.0));
			const Vector k3 = f(t + (h * 2.0 / 3.0), x + h * (-k1 / 3.0 + k2));
			const Vector k4 = f(t + h, x + h * (k1 - k2 + k1));

			return x + (k1 + 3.0 * k2 + 3.0 * k3 + k4) * (h / 8.0);
		}


		/// Compute one step of the Adams-Bashforth linear multistep method of
		/// 2nd order for ordinary differential equations.
		/// This function is used in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_adams2(
			OdeFunction f, const Vector& x0, real t0,
			const Vector& x1, real t1, real h = 0.001) {

			return x1 + h * ((3. / 2.) * f(t1, x1) - f(t0, x0) / 2.);
		}


		/// Compute one step of the Adams-Bashforth linear multistep method of
		/// 3rd order for ordinary differential equations.
		/// This function is used in solvers to solve an ODE over a certain domain.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x The starting vector of variables
		/// @param t The starting value of the time (independent variable)
		/// @param h The step size
		/// @return The resulting vector of variables
		template<typename Vector, typename OdeFunction = ode_function<Vector>>
		inline Vector step_adams3(
			OdeFunction f, const Vector& x0, real t0, const Vector& x1, real t1,
			const Vector& x2, real t2, real h = 0.001) {

			return x2 + h *	(
				(23. / 12.) * f(t2, x2) - (4. / 3.) * f(t1, x1) + (3. / 8.) * f(t0, x0)
			);
		}


		// Solvers
		// (functions which solve numerically an ODE over an interval)


		/// Integrate an ordinary differential equation using any numerical algorithm
		/// with a constant step size, such as Runge-Kutta methods. This function
		/// does not use a specific method but uses the step argument function
		/// to iterate each step of an arbitrary fixed step algorithm. If the step size
		/// does not exactly cover the interval of integration, the last step is shortened.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param step A function which integrates numerically the differential equation
		/// between \f$t\f$ and \f$t + h\f$, such as the functions named ode::step_*
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>,
			typename StepFunction = std::function<Vector(OdeFunction, real, const Vector&)>
		>
		inline ode_solution_t<Vector>
		solve_fixstep(
			OdeFunction f, const Vector& x0, real t0, real tf,
			StepFunction step, real stepsize = 0.001) {

			if (tf < t0) {
				TH_MATH_ERROR("ode::solve_fixstep", tf, INVALID_ARGUMENT);
				ode_solution_t<Vector> err; err.t = Vector(nan());
				return err;
			}

			const unsigned int steps = floor((tf - t0) / stepsize);
			unsigned int total_steps = steps;
			unsigned int i;

			if (abs(t0 + steps * stepsize - tf) > MACH_EPSILON)
				total_steps++;

			// Initialize solution structure
			auto solution = ode_solution_t<Vector>(total_steps + 1, x0, t0);


			// Iterate over each step of the numerical method
			for (i = 1; i <= steps; ++i) {
				solution.x[i] = step(f, solution.x[i - 1], solution.t[i - 1], stepsize);
				solution.t[i] = solution.t[i - 1] + stepsize;
			}


			// Additional shorter step if the stepsize
			// does not cover exactly the time interval
			if (steps != total_steps) {
				solution.x[i] = step(
					f, solution.x[i - 1], solution.t[i - 1], tf - solution.t[i - 1]);
				solution.t[i] = tf;
			}

			return solution;
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using Euler's method.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_euler(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.0001) {

			return solve_fixstep(f, x0, t0, tf, step_euler<Vector, OdeFunction>, stepsize);
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using the midpoint method.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_midpoint(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.0001) {

			return solve_fixstep(f, x0, t0, tf, step_midpoint<Vector, OdeFunction>, stepsize);
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using Heun's method.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_heun(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.0001) {

			return solve_fixstep(f, x0, t0, tf, step_heun<Vector, OdeFunction>, stepsize);
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using Runge-Kutta's method
		/// of 2nd order.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_rk2(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.0001) {

			return solve_fixstep(f, x0, t0, tf, step_rk2<Vector, OdeFunction>, stepsize);
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using Runge-Kutta's method
		/// of 4th order.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_rk4(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.01) {

			return solve_fixstep(f, x0, t0, tf, step_rk4<Vector, OdeFunction>, stepsize);
		}


		/// Integrate an ordinary differential equation over a certain
		/// domain with the given initial conditions using Kutta's 3/8 rule method.
		///
		/// @param f A function representing the system of differential equations,
		/// following the signature of ode_function.
		/// @param x0 The initial value of the variables
		/// @param t0 The starting value of the time variable
		/// @param tf The final value of the time variable
		/// @param stepsize The constant step size
		/// @return The numerical solution of the equation, as an ode_solution_t
		/// structure, holding a vector t of time values and a vector x of the variables.
		template <
			typename Vector, typename OdeFunction = ode_function<Vector>
		>
		inline ode_solution_t<Vector> solve_k38(
			OdeFunction f, const Vector& x0,
			real t0, real tf, real stepsize = 0.0001) {

			return solve_fixstep(f, x0, t0, tf, step_k38<Vector, OdeFunction>, stepsize);
		}
	}
}

#endif
