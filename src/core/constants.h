
///
/// @file constants.h Mathematical constants and default algorithm
/// parameters. You may change the library's default behavior
/// by redefining the macros starting with THEORETICA_
///

#ifndef THEORETICA_CONSTANTS_H
#define THEORETICA_CONSTANTS_H


#include <limits>
#include <cstdint>
#include <type_traits>


/// THEORETICA_DISABLE_X86 Define this macro to disable
/// Assembly x86 optimizations.
#ifndef THEORETICA_DISABLE_X86

#ifndef THEORETICA_X86
#if defined(__amd64__) || defined(__amd64) || defined(__x86_64__) \
	|| defined(_M_AMD64) || defined(_M_X64) || defined(i386) \
	|| defined(__i386) || defined(__i386__) || defined(_M_IX86) \
	|| defined(_X86_) || defined(_M_I86) || defined(__X86__)

/// THEORETICA_X86 This macro is automatically defined
/// on most compilers to enable Assembly x86 optimizations.
/// @see THEORETICA_DISABLE_X86
#define THEORETICA_X86
#endif
#endif
#endif


/// Order of Taylor series approximations
#ifndef THEORETICA_TAYLOR_ORDER
#define THEORETICA_TAYLOR_ORDER 12
#endif

/// Relative precision for derivative approximation
#ifndef THEORETICA_DERIV_PREC
#define THEORETICA_DERIV_PREC 10E+5
#endif

/// Default number of steps for integral approximation
#ifndef THEORETICA_INTEGRATION_STEPS
#define THEORETICA_INTEGRATION_STEPS 100
#endif

#ifndef THEORETICA_RKDP_TOL
#define THEORETICA_RKDP_TOL 1E-8
#endif

// Default tolerance for integral approximation
#ifndef THEORETICA_INTEGRATION_TOL
#define THEORETICA_INTEGRATION_TOL 1E-08
#endif

/// Biggest fractional part to ignore in powf computation
#ifndef THEORETICA_POWF_APPROX_TOL
#define THEORETICA_POWF_APPROX_TOL 1E-08
#endif
	
/// Approximation tolerance for root finding
#ifndef THEORETICA_ROOT_APPROX_TOL
#define THEORETICA_ROOT_APPROX_TOL 1E-08
#endif

/// Approximation tolerance for Newton's method root finding
#ifndef THEORETICA_NEWTON_RAPHSON_TOL
#define THEORETICA_NEWTON_RAPHSON_TOL 1E-08
#endif

/// Approximation tolerance for bisection root finding
#ifndef THEORETICA_BISECTION_APPROX_TOL
#define THEORETICA_BISECTION_APPROX_TOL 1E-08
#endif

/// Maximum number of iterations for bisection
#ifndef THEORETICA_MAX_BISECTION_ITER
#define THEORETICA_MAX_BISECTION_ITER 100
#endif

/// Maximum number of iterations for golden section search
#ifndef THEORETICA_MAX_GOLDENSECTION_ITER
#define THEORETICA_MAX_GOLDENSECTION_ITER 100
#endif

/// Maximum number of iterations for Halley's method
#ifndef THEORETICA_MAX_HALLEY_ITER
#define THEORETICA_MAX_HALLEY_ITER 100
#endif

/// Maximum number of iterations for Newton-Raphson root finding
#ifndef THEORETICA_MAX_NEWTON_ITER
#define THEORETICA_MAX_NEWTON_ITER 100
#endif

/// Maximum number of iterations for Steffensen root finding
#ifndef THEORETICA_MAX_STEFFENSEN_ITER
#define THEORETICA_MAX_STEFFENSEN_ITER 100
#endif

/// Maximum number of iterations for Chebyshev root finding
#ifndef THEORETICA_MAX_CHEBYSHEV_ITER
#define THEORETICA_MAX_CHEBYSHEV_ITER 100
#endif

/// Maximum number of failed iterations for the Try-and-Catch algorithm
#ifndef THEORETICA_MAX_TRYANDCATCH_ITER
#define THEORETICA_MAX_TRYANDCATCH_ITER 100
#endif

/// Default variation for derivative approximation
#ifndef THEORETICA_DERIV_STEPSIZE
#define THEORETICA_DERIV_STEPSIZE 1E-3
#endif

/// Default step size for gradient descent minimization
#ifndef THEORETICA_MINGRAD_GAMMA
#define THEORETICA_MINGRAD_GAMMA -0.005
#endif

/// Default tolerance for gradient descent minimization
#ifndef THEORETICA_MINGRAD_TOLERANCE
#define THEORETICA_MINGRAD_TOLERANCE 1E-3
#endif

/// Maximum number of iterations for gradient descent minimization
#ifndef THEORETICA_MINGRAD_MAX_ITER
#define THEORETICA_MINGRAD_MAX_ITER 50000
#endif


/// Maximum number of polynomial division iterations
#ifndef THEORETICA_MAX_POLYNDIV_ITER
#define THEORETICA_MAX_POLYNDIV_ITER 100
#endif


#ifndef THEORETICA_RAND_REAL_PREC

/// Default precision for random number generation using rand_uniform()
#ifdef THEORETICA_FLOAT_PREC
#define THEORETICA_RAND_REAL_PREC (uint64_t(1) << 23)
#else
#define THEORETICA_RAND_REAL_PREC (uint64_t(1) << 31)
#endif

#endif


/// Default depth of the Metropolis algorithm
#ifndef THEORETICA_METROPOLIS_DEPTH
#define THEORETICA_METROPOLIS_DEPTH 16
#endif


/// Enable constexpr in function declarations if C++14 is supported.
#if (__cplusplus >= 201402L)
#define TH_CONSTEXPR constexpr
#else
#define TH_CONSTEXPR
#endif


/// Enable constexpr in if statements if C++17 is supported.
#if (__cplusplus >= 201703L)
#define TH_CONSTIF constexpr
#else
#define TH_CONSTIF
#endif


/// @namespace theoretica Main namespace of the library which contains all functions and objects
namespace theoretica {


	/// A real number, defined as a floating point type.
	///
	/// The underlying type is determined by the defined macros:
	/// By default, `real` will be defined as the `double` type.
	/// If `THEORETICA_FLOAT_PREC` is defined, `real` will be defined as a `float`,
	/// if `THEORETICA_LONG_DOUBLE_PREC` is defined, `real` will be defined as a `long double`
	/// @note The `THEORETICA_ARBITRARY_PREC` option is currently unsupported

#ifdef THEORETICA_LONG_DOUBLE_PREC

	using real = long double;

#elif defined(THEORETICA_FLOAT_PREC)

	using real = float;

#elif defined(THEORETICA_ARBITRARY_PREC)

// TO-DO bigfloat arbitrary precision

#else

	using real = double;

#endif


	// Type traits.


	// Implement a simple void_t trait in C++14.
	namespace _internal {
		template<typename ...Args>
		struct make_void { typedef void type; };
		template<typename ...Args>
		using void_t = typename make_void<Args...>::type;
	}


	/// Type trait to check whether a type represents
	/// a real number.
	template<typename Type>
	struct is_real_type : std::is_floating_point<Type> {};


	/// Type trait to check whether a type represents
	/// a real number.
	template<>
	struct is_real_type<real> : std::true_type {};


	/// Extract the type of an indexable container from its operator[].
	template<typename Structure>
	using get_indexable_element_t =
		std::remove_reference_t<decltype(std::declval<Structure>()[0])>;


	/// Type trait to check whether an indexable container
	/// has elements of the given type.
	template<typename Structure, typename Type>
	struct has_type_elements
	: std::is_same<get_indexable_element_t<Structure>, Type> {};


	/// Type trait to check whether an indexable container
	/// has real elements.
	template<typename Structure>
	using has_real_elements = is_real_type<get_indexable_element_t<Structure>>;


	/// Check whether a structure is orderable,
	/// by checking that it has a comparison operator<().
	template<typename Structure, typename = _internal::void_t<>>
	struct is_orderable : std::false_type{};


	/// Check whether a structure is orderable,
	/// by checking that it has a comparison operator<().
	template<typename Structure>
	struct is_orderable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>() < std::declval<Structure>())>
	> : std::true_type{};


	/// Check whether a structure is indexable by a single integer index,
	/// by checking that it has the operator[](0).
	template<typename Structure, typename = _internal::void_t<>>
	struct is_indexable : std::false_type{};


	/// Check whether a structure is indexable by a single integer index,
	/// by checking that it has the operator[](0).
	template<typename Structure>
	struct is_indexable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>()[0])>
	> : std::true_type{};


	/// Check whether a structure is iterable,
	/// by checking that it has a method begin().
	template<typename Structure, typename = _internal::void_t<>>
	struct is_iterable : std::false_type{};


	/// Check whether a structure is iterable,
	/// by checking that it has a method begin().
	template<typename Structure>
	struct is_iterable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>().begin())>
	> : std::true_type{};



	// Mathematical constants and default algorithm parameters.


	/// Machine epsilon for the real type
	constexpr real MACH_EPSILON = std::numeric_limits<real>::epsilon();

	/// The Phi (Golden Section) mathematical constant
	constexpr real PHI = 1.6180339887498948482045868;

	/// The inverse of the Golden Section mathematical constant
	constexpr real INVPHI = 0.6180339887498948482045868;

	/// The Pi mathematical constant
	constexpr real PI = 3.141592653589793238462643;

	/// Half of Pi
	constexpr real PI2 = 1.57079632679489655799898;

	/// A quarter of Pi
	constexpr real PI4 = PI / 4.0;

	/// Pi multiplied by 2
	constexpr real PIDOUBLE = PI * 2;

	/// The Tau mathematical constant (Pi times 2)
	constexpr real TAU = PI * 2;

	/// The inverse of Pi
	constexpr real INVPI = 1.0 / PI;

	/// The square root of Pi
	constexpr real SQRTPI = 1.7724538509055159927;

	/// The Euler mathematical constant (e)
	constexpr real E = 2.718281828459045235360287;

	/// The binary logarithm of e
	constexpr real LOG2E = 1.44269504088896338700465094;

	/// The binary logarithm of 10
	constexpr real LOG210 = 3.32192809488736218170856773213;

	/// The base-10 logarithm of e
	constexpr real LOG10E = 0.434294481903;

	/// The natural logarithm of 2
	constexpr real LN2 = 0.69314718056;

	/// The natural logarithm of 10
	constexpr real LN10 = 2.30258509299;

	/// The scalar conversion factor from degrees to radians
	constexpr real DEG2RAD = 0.017453292519943295474371680598;

	/// The scalar conversion factor from radians to degrees
	constexpr real RAD2DEG = 57.2957795130823228646477218717;

	/// The square root of 2
	constexpr real SQRT2 = 1.4142135623730950488;

	/// The inverse of the square root of 2
	constexpr real INVSQR2 = 0.7071067811865475;

	/// The square root of 3
	constexpr real SQRT3 = 1.732050807568877;

	/// Order of Taylor series approximations
	constexpr int TAYLOR_ORDER = THEORETICA_TAYLOR_ORDER;

	/// Default number of steps for integral approximation
	constexpr int INTEGRATION_STEPS = THEORETICA_INTEGRATION_STEPS;

	// Default tolerance for integral approximation
	constexpr real INTEGRATION_TOL = THEORETICA_INTEGRATION_TOL;

	/// Relative precision for derivative approximation
	constexpr real DERIV_PREC = THEORETICA_DERIV_PREC;

	/// Biggest fractional part to ignore in powf computation
	constexpr real POWF_POWER_TOLERANCE = THEORETICA_POWF_APPROX_TOL;

	/// Approximation tolerance for root finding
	constexpr real ROOT_APPROX_TOL = THEORETICA_ROOT_APPROX_TOL;

	/// Approximation tolerance for the bisection algorithm
	constexpr real BISECTION_APPROX_TOL = THEORETICA_BISECTION_APPROX_TOL;

	/// Approximation tolerance for the Newton-Raphson algorithm
	constexpr real NEWTON_RAPHSON_TOL = THEORETICA_NEWTON_RAPHSON_TOL;

	/// Maximum number of iterations for the bisection algorithm
	constexpr unsigned int MAX_BISECTION_ITER = THEORETICA_MAX_BISECTION_ITER;

	/// Maximum number of iterations for the golden section search algorithm
	constexpr unsigned int MAX_GOLDENSECTION_ITER = THEORETICA_MAX_GOLDENSECTION_ITER;

	/// Maximum number of iterations for Halley's method
	constexpr unsigned int MAX_HALLEY_ITER = THEORETICA_MAX_HALLEY_ITER;

	/// Maximum number of iterations for the Newton-Raphson algorithm
	constexpr unsigned int MAX_NEWTON_ITER = THEORETICA_MAX_NEWTON_ITER;

	/// Maximum number of iterations for the Steffensen algorithm
	constexpr unsigned int MAX_STEFFENSEN_ITER = THEORETICA_MAX_STEFFENSEN_ITER;

	/// Maximum number of iterations for the Chebyshev algorithm
	constexpr unsigned int MAX_CHEBYSHEV_ITER = THEORETICA_MAX_CHEBYSHEV_ITER;

	/// Maximum number of failed iterations for the Try-and-Catch algorithm
	constexpr unsigned int MAX_TRYANDCATCH_ITER = THEORETICA_MAX_TRYANDCATCH_ITER;

	/// Default variation for derivative approximation
	constexpr real DERIV_STEPSIZE = THEORETICA_DERIV_STEPSIZE;

	/// Default step size for gradient descent minimization
	constexpr real MINGRAD_GAMMA = THEORETICA_MINGRAD_GAMMA;

	/// Default tolerance for gradient descent minimization
	constexpr real MINGRAD_TOLERANCE = THEORETICA_MINGRAD_TOLERANCE;

	/// Maximum number of iterations for gradient descent minimization
	constexpr unsigned int MINGRAD_MAX_ITER = THEORETICA_MINGRAD_MAX_ITER;

	/// Default precision for random number generation using rand_uniform()
	constexpr uint64_t RAND_REAL_PREC = THEORETICA_RAND_REAL_PREC;

	/// Default depth of the Metropolis algorithm
	constexpr unsigned int METROPOLIS_DEPTH = THEORETICA_METROPOLIS_DEPTH;

}

#endif
