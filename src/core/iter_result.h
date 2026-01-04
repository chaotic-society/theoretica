
///
/// @file iter_result.h Structured results for iterative algorithms, containing information about convergence.
///

#ifndef THEORETICA_ITER_RESULT_H
#define THEORETICA_ITER_RESULT_H

#include <string>
#include "./constants.h"
#include "./error.h"


namespace theoretica {


	/// @enum ConvergenceStatus
	/// Status codes for iterative algorithm termination
	enum class ConvergenceStatus {

		/// Algorithm converged successfully
		Success,

		/// Maximum iterations exceeded
		MaxIterations,

		/// No progress in iterations
		Stalled,

		/// Invalid input provided
		InvalidInput,   
		
		/// Algorithm diverged
		Diverged,

		/// User terminated early
		UserInterrupt
	};


	/// @class iter_result
	/// A structure returned by iterative algorithms containing
	/// the computed value, convergence information, and diagnostics.
	///
	/// Features implicit conversion to the result type.
	/// @tparam T The type of the computed result (usually real or vec<real, N>)
	template<typename Type = real>
	struct iter_result {

		/// The computed result value
		Type value;

		/// Whether the algorithm converged within the specified criteria
		bool converged = false;

		/// Status code indicating reason for termination
		ConvergenceStatus status = ConvergenceStatus::Success;

		/// Number of iterations performed
		unsigned int iterations = 0;

		/// Final error or residual norm (exact meaning depends on the algorithm)
		real residual = inf();


		/// Constructor with default values
		iter_result() : value(Type(nan())) {}


		/// Constructor with value for reporting success
		iter_result(const Type& val) : value(val) {
			converged = true;
			status = ConvergenceStatus::Success;
		}
		
		iter_result(const Type& value, unsigned int iterations)
		: value(value), iterations(iterations) {
			converged = true;
			status = ConvergenceStatus::Success;
		}

		iter_result(const Type& value, unsigned int iterations, real residual)
		: value(value), iterations(iterations), residual(residual) {
			converged = true;
			status = ConvergenceStatus::Success;
		}


		/// Constructor with convergence status for reporting failure
		iter_result(ConvergenceStatus status, unsigned int iterations = 0)
		: status(status), iterations(iterations) {
			value = Type(nan());
			converged = false;
		}

		iter_result(ConvergenceStatus status, unsigned int iterations, real residual)
		: status(status), iterations(iterations), residual(residual) {
			value = Type(nan());
			converged = false;
		}


		/// Implicit conversion to result type, allows code like: Type res = algorithm();
		operator Type() const {
			return value;
		}

		/// Explicit conversion to boolean (converged status), allows code like: if (result) { ... }
		explicit operator bool() const {
			return converged;
		}


		/// Get a human-readable string description of the status of convergence of an iterative algorithm.
		inline std::string status_string() const {

			switch(status) {
				case ConvergenceStatus::Success: 
					return "Converged successfully";
				case ConvergenceStatus::MaxIterations:
					return "Maximum iterations exceeded without converging to desired accuracy";
				case ConvergenceStatus:: Stalled:
					return "Algorithm stalled";
				case ConvergenceStatus::InvalidInput:
					return "Invalid input provided";
				case ConvergenceStatus::Diverged:
					return "Algorithm diverged";
				case ConvergenceStatus::UserInterrupt: 
					return "User interrupt";
				default:
					return "Unknown status";
			}
		}

		
#ifndef THEORETICA_NO_PRINT

		/// Convert the iter_result number to string representation
		inline std::string to_string() const {

			std::stringstream res;

			res << "Value = " << value << "\n";

			if (!converged)
				res << "Converged = " << (converged ? "true" : "false") << "\n";

			res << "Status: " << status_string() << "\n";
			res << "Iterations = " << iterations << "\n";
			res << "Residual = " << residual;

			return res.str();
		}


		/// Stream the complex number in string representation
		/// to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const iter_result<Type>& obj) {
			return out << obj.to_string();
		}

#endif

	};

}

#endif
