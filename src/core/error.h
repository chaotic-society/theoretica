
///
/// @file core/error.h Error handling
///

#ifndef THEORETICA_ERROR_H
#define THEORETICA_ERROR_H

#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)
#include <exception>
#include <string>

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#endif

#endif

#include <cerrno>
#include <limits>
#include "./constants.h"


namespace theoretica {

	/// Math error enumeration
	enum MATH_ERRCODE{
		NO_ERROR = 0x00, // No error
		DIV_BY_ZERO = 0x01, // Division by zero
		OUT_OF_DOMAIN = 0x02, // An argument is out of range
		OUT_OF_RANGE = 0x04, // The result would be out of range
		IMPOSSIBLE_OPERATION = 0x08, // Impossible operation
		NO_ALGO_CONVERGENCE = 0x10, // The algorithm did not converge
		INVALID_ARGUMENT = 0x20 // Invalid argument size or value
	};


	/// Convert a MATH_ERRCODE to errno error codes
	inline int th_errcode_to_errno(MATH_ERRCODE err) {
		switch(err) {
			case NO_ERROR: return 0; break;
			case DIV_BY_ZERO: return ERANGE; break;
			case OUT_OF_DOMAIN: return EDOM; break;
			case OUT_OF_RANGE: return ERANGE; break;
			case IMPOSSIBLE_OPERATION: return EDOM; break;
			case NO_ALGO_CONVERGENCE: return ERANGE; break;
			case INVALID_ARGUMENT: return EINVAL; break;
			default: return 0; break;
		}
	}


	/// Return a quiet NaN number in floating point representation
	inline real nan() {
		return std::numeric_limits<real>::quiet_NaN();
	}


	/// Check whether a generic variable is
	/// (equivalent to) a NaN number.
	///
	/// NaN numbers are the only variables which do not
	/// compare equal to themselves in floating point
	/// operations. This is valid for real types but
	/// also for any mathematical structure, as NaNs
	/// are used to report failure inside the library.
	template<typename T>
	inline bool is_nan(const T& x) {
		return !(x == x);
	}


	/// Return positive infinity in floating point representation
	inline real inf() {
		return std::numeric_limits<real>::infinity();
	}


#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)

	class math_exception : std::exception {

	private:
		MATH_ERRCODE err;
		std::string func_name;
		std::string file_name;
		unsigned int code_line;
		real val;

	public:
		math_exception(MATH_ERRCODE a_err, const std::string& a_func_name,
			const std::string& a_file_name, unsigned int a_code_line, real a_val)
				: err(a_err), func_name(a_func_name), file_name(a_file_name),
					code_line(a_code_line), val(a_val) {}

		~math_exception() = default;


		/// Return a string describing the exception
		inline const char* what() const noexcept {

			switch(err) {
				case NO_ERROR: return "No error"; break;
				case DIV_BY_ZERO: return "Division by zero"; break;
				case OUT_OF_DOMAIN:
					return "An argument was out of the domain of the called function"; break;
				case IMPOSSIBLE_OPERATION:
					return "A mathematically impossible operation was requested"; break;
				case NO_ALGO_CONVERGENCE: return "The algorithm did not converge"; break;
				case INVALID_ARGUMENT: return "Invalid argument size or value"; break;
				default: return "Unknown error"; break;
			}
		}


		/// Get the error code associated with the exception
		inline MATH_ERRCODE err_code() const {
			return err;
		}


		/// Get the name of the throwing function
		inline std::string get_function_name() const {
			return func_name;
		}


		/// Get the name of the file in which the exception was thrown
		inline std::string get_file_name() const {
			return file_name;
		}


		/// Get the line number at which the exception was thrown
		inline unsigned int get_line_number() const {
			return code_line;
		}


		/// Get a real value associated with the exception
		inline real get_value() const {
			return val;
		}
		

#ifndef THEORETICA_NO_PRINT

		/// Get a string representation of the exception
		inline std::string to_string() const {

			std::stringstream err_str;

			err_str << file_name << "(" << code_line << "):";
			err_str << func_name << "(" << val << "): ";

			switch(err) {
				case NO_ERROR: err_str << "No error"; break;
				case DIV_BY_ZERO: err_str << "Division by zero"; break;
				case OUT_OF_DOMAIN:
					err_str << "An argument was out of the domain of the called function"; break;
				case IMPOSSIBLE_OPERATION:
					err_str << "A mathematically impossible operation was requested"; break;
				case NO_ALGO_CONVERGENCE: err_str << "The algorithm did not converge"; break;
				case INVALID_ARGUMENT: err_str << "Invalid argument size or value"; break;
				default: err_str << "Unknown error"; break;
			}

			return err_str.str();
		}


		/// Stream the exception in string representation
		/// to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const math_exception& obj) {
			return out << obj.to_string();
		}

#endif

	};

#endif

}


/// TH_MATH_ERROR is a macro which throws exceptions
/// or modifies errno (depending on which compiling options
/// are defined)

	
// Only throw exceptions, without modifying errno
#ifdef THEORETICA_ONLY_EXCEPTIONS

#define TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	throw math_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE);

#define TH_MATH_ERROR_R(F_NAME, VALUE, EXCEPTION) \
	TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION)

// Throw exceptions and modify errno
#elif defined(THEORETICA_THROW_EXCEPTIONS)

#define TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	errno = th_errcode_to_errno(EXCEPTION); \
	throw math_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE);

// Modify errno only by default
#else

#define TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	errno = th_errcode_to_errno(EXCEPTION);

#endif


#endif
