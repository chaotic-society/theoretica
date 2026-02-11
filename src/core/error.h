
///
/// @file core/error.h Error handling
///

#ifndef THEORETICA_ERROR_H
#define THEORETICA_ERROR_H

#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)

#include <exception>

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#endif

#endif

#include <string>
#include <cerrno>
#include <limits>
#include "./constants.h"


namespace theoretica {


	/// Math error enumeration
    enum class MathError : int {
        None = 0x00,				///< No error
        DivByZero = 0x01,			///< Division by zero
        OutOfDomain = 0x02,			///< Argument out of domain
        OutOfRange = 0x04,			///< Result out of range
        ImpossibleOperation = 0x08,	///< Mathematically impossible operation
        NoConvergence = 0x10,		///< Algorithm did not converge
        InvalidArgument = 0x20		///< Invalid argument
    };


	/// Convert an MathError class enum to conventional errno codes.
	inline int to_errno(MathError err) {
		switch(err) {
			case MathError::None: return 0; break;
			case MathError::DivByZero: return ERANGE; break;
			case MathError::OutOfDomain: return EDOM; break;
			case MathError::OutOfRange: return ERANGE; break;
			case MathError::ImpossibleOperation: return EDOM; break;
			case MathError::NoConvergence: return ERANGE; break;
			case MathError::InvalidArgument: return EINVAL; break;
			default: return 0; break;
		}
	}


	/// Convert a MathError class enum to a string description.
	inline const char* to_cstring(MathError err) {
		switch (err) {
			case MathError::None: return "No error"; break;
			case MathError::DivByZero: return "Division by zero"; break;
			case MathError::OutOfDomain:
				return "An argument was out of the domain of the called function"; break;
			case MathError::ImpossibleOperation:
				return "A mathematically impossible operation was requested"; break;
			case MathError::NoConvergence: return "The algorithm did not converge"; break;
			case MathError::InvalidArgument: return "Invalid argument size or value"; break;
			default: return "Unknown error"; break;
		}
	}


	/// Convert a MathError class enum to a string description.
	inline std::string to_string(MathError err) {
		return to_cstring(err);
	}


	/// Return a quiet NaN number in floating point representation
	inline TH_CONSTEXPR real nan() {
		return std::numeric_limits<real>::quiet_NaN();
	}


	/// Check whether a generic variable is (equivalent to) a NaN number.
	///
	/// NaN numbers are the only variables which do not
	/// compare equal to themselves in floating point
	/// operations. This is valid for real types but
	/// also for any mathematical structure, as NaNs
	/// are used to report failure inside the library.
	///
	/// @param x The mathematical structure to test
	/// for being a NaN or NaN-equivalent structure.
	template<typename T>
	inline bool is_nan(const T& x) {
		return !(x == x);
	}


	/// Get positive infinity in floating point representation
	inline TH_CONSTEXPR real inf() {
		return std::numeric_limits<real>::infinity();
	}


	/// Check whether a real number is infinite.
	inline bool is_inf(real x) {
		return (x == inf()) || (x == -inf());
	}


#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)

	class math_exception : std::exception {

	private:
		MathError err;
		std::string func_name;
		std::string file_name;
		unsigned int code_line;
		real val;

	public:
		math_exception(MathError a_err, const std::string& a_func_name,
			const std::string& a_file_name, unsigned int a_code_line, real a_val)
				: err(a_err), func_name(a_func_name), file_name(a_file_name),
					code_line(a_code_line), val(a_val) {}

		~math_exception() = default;


		/// Return a string describing the exception
		inline const char* what() const noexcept {
			return to_cstring(err);
		}


		/// Get the error code associated with the exception
		inline MathError err_code() const {
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
				case MathError::None: err_str << "No error"; break;
				case MathError::DivByZero: err_str << "Division by zero"; break;
				case MathError::OutOfDomain:
					err_str << "An argument was out of the domain of the called function"; break;
				case MathError::ImpossibleOperation:
					err_str << "A mathematically impossible operation was requested"; break;
				case MathError::NoConvergence: err_str << "The algorithm did not converge"; break;
				case MathError::InvalidArgument: err_str << "Invalid argument size or value"; break;
				default: err_str << "Unknown error"; break;
			}

			return err_str.str();
		}


		/// Convert the exception to string representation.
		inline operator std::string() {
			return to_string();
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
	{ throw theoretica::math_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE); }

#define TH_MATH_ERROR_R(F_NAME, VALUE, EXCEPTION) \
	TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION)

// Throw exceptions and modify errno
#elif defined(THEORETICA_THROW_EXCEPTIONS)

#define TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	{ errno = theoretica::to_errno(EXCEPTION); \
	throw math_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE); }

// Modify errno only by default
#else

#define TH_MATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	{ errno = theoretica::to_errno(EXCEPTION); }

#endif


// Output the value of an expression with additional information,
// for debugging purposes.
#define TH_DEBUG(VARIABLE) { \
	std::cout << __FILE__ << ":" << __LINE__ << ": " \
	<< #VARIABLE << " = " << VARIABLE << std::endl; }


#endif
