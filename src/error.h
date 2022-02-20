#ifndef UROBORO_ERROR_H
#define UROBORO_ERROR_H

#if defined(UROBORO_THROW_EXCEPTIONS) || defined(UROBORO_ONLY_EXCEPTIONS)
#include <exception>
#include <string>
#include <sstream>
#endif

#include <cerrno>
#include <limits>
#include "./constants.h"


namespace uroboro {

	// Math error enumeration
	enum UMATH_ERRCODE{
		NO_ERROR = 0x00, // No error
		DIV_BY_ZERO = 0x01, // Division by zero
		OUT_OF_DOMAIN = 0x02, // An argument is out of range
		OUT_OF_RANGE = 0x04, // The result would be out of range
		IMPOSSIBLE_OPERATION = 0x08, // Impossible operation
		NO_ALGO_CONVERGENCE = 0x10, // The algorithm did not converge
		INVALID_ARGUMENT = 0x20 // Invalid argument size or value
	};


	// Convert a UMATH_ERRCODE to errno error codes
	int umath_errcode_to_errno(UMATH_ERRCODE err) {
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


	// Return a quiet NaN number
	inline real nan() {
		return std::numeric_limits<real>::quiet_NaN();
	}


	// Return positive infinity
	inline real inf() {
		return std::numeric_limits<real>::infinity();
	}


#if defined(UROBORO_THROW_EXCEPTIONS) || defined(UROBORO_ONLY_EXCEPTIONS)

	class MathException : std::exception {

	private:
		UMATH_ERRCODE err;
		std::string func_name;
		std::string file_name;
		unsigned int code_line;
		real val;

	public:
		MathException(UMATH_ERRCODE a_err, std::string a_func_name,
			std::string a_file_name, unsigned int a_code_line, real a_val) {

			err = a_err;
			func_name = a_func_name;
			file_name = a_file_name;
			code_line = a_code_line;
			val = a_val;
		}

		~MathException() = default;


		// Return a string describing the exception
		inline const char* what() {

			// std::string w = std::string(file_name);

			// w += ":";
			// w += func_name;
			// w += ": ";

			switch(err) {
				case UMATH_ERRCODE::NO_ERROR: return "No error"; break;
				case UMATH_ERRCODE::DIV_BY_ZERO: return "Division by zero"; break;
				case UMATH_ERRCODE::OUT_OF_DOMAIN:
					return "An argument was out of the domain of the called function"; break;
				case UMATH_ERRCODE::IMPOSSIBLE_OPERATION:
					return "A mathematically impossible operation was requested"; break;
				case UMATH_ERRCODE::NO_ALGO_CONVERGENCE: return "The algorithm did not converge"; break;
				case UMATH_ERRCODE::INVALID_ARGUMENT: return "Invalid argument size or value"; break;
				default: return "Unknown error"; break;
			}

			// return w.c_str();
		}


		inline UMATH_ERRCODE err_code() {
			return err;
		}

	};

#endif

}


// UMATH_ERROR is a macro which throws exceptions
// or modifies errno depending on which compiling options
// are defined

	
// Only throw exceptions, without modifying errno
#ifdef UROBORO_ONLY_EXCEPTIONS

#define UMATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	throw MathException(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE);

#define UMATH_ERROR_R(F_NAME, VALUE, EXCEPTION) \
	UMATH_ERROR(F_NAME, VALUE, EXCEPTION)

// Throw exceptions and modify errno
#elif defined(UROBORO_THROW_EXCEPTIONS)

#define UMATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	errno = umath_errcode_to_errno(EXCEPTION); \
	throw MathException(EXCEPTION, F_NAME, __FILE__, __LINE__, VALUE);

#define UMATH_ERROR_R(F_NAME, VALUE, EXCEPTION) \
	UMATH_ERROR(F_NAME, VALUE, EXCEPTION)

// Modify errno only
#else

// For a generic function (NO return statement)
#define UMATH_ERROR(F_NAME, VALUE, EXCEPTION) \
	errno = umath_errcode_to_errno(EXCEPTION);


// For a real function
#define UMATH_ERROR_R(F_NAME, VALUE, EXCEPTION) \
	errno = umath_errcode_to_errno(EXCEPTION); \
	return uroboro::nan();

#endif


#endif
