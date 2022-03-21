#ifndef UROBORO_CONSTANTS_H
#define UROBORO_CONSTANTS_H


// Order of Taylor series approximations
#ifndef UROBORO_TAYLOR_ORDER
#define UROBORO_TAYLOR_ORDER 7
#endif

// Relative precision for derivative approximation
#ifndef UROBORO_DERIV_PREC
#define UROBORO_DERIV_PREC 100000.0
#endif

// Default number of steps for integral approximation
#ifndef UROBORO_INTEGRATION_STEPS
#define UROBORO_INTEGRATION_STEPS 100
#endif

// Biggest fractional part to ignore in powf computation
#ifndef UROBORO_POWF_APPROX_TOL
#define UROBORO_POWF_APPROX_TOL 0.00000001
#endif
	
// Approximation tolerance for root finding
#ifndef UROBORO_ROOT_APPROX_TOL
#define UROBORO_ROOT_APPROX_TOL 0.00000001
#endif


// Approximation tolerance for root finding
#ifndef UROBORO_NEWTON_RAPHSON_TOL
#define UROBORO_NEWTON_RAPHSON_TOL 0.00000001
#endif

// Approximation tolerance for bisection root finding
#ifndef UROBORO_BISECTION_APPROX_TOL
#define UROBORO_BISECTION_APPROX_TOL 0.00000001
#endif

// Maximum number of iterations for bisection
#ifndef UROBORO_MAX_BISECTION_ITER
#define UROBORO_MAX_BISECTION_ITER 100
#endif

// Maximum number of iterations for Newton-Raphson
#ifndef UROBORO_MAX_NEWTON_ITER
#define UROBORO_MAX_NEWTON_ITER 100
#endif

// Maximum number of iterations for Steffensen
#ifndef UROBORO_MAX_STEFFENSEN_ITER
#define UROBORO_MAX_STEFFENSEN_ITER 100
#endif

// Maximum number of iterations for Chebyshev
#ifndef UROBORO_MAX_CHEBYSHEV_ITER
#define UROBORO_MAX_CHEBYSHEV_ITER 100
#endif


namespace uroboro {


	// Real number type

#ifdef UROBORO_LONG_DOUBLE_PREC

	using real = long double;

#elif defined(UROBORO_FLOAT_PREC)

	using real = float;

#elif defined(UROBORO_ARBITRARY_PREC)

// TO-DO bigfloat arbitrary precision

#else

	using real = double;

#endif


	// PI constant.
	constexpr real PI = 3.141592653589793238462643;

	// Half PI.
	constexpr real PI2 = PI / 2.0;

	// A quarter of PI.
	constexpr real PI4 = PI / 4.0;

	// PI multiplied by 2.
	constexpr real PIDOUBLE = PI * 2;

	// TAU constant (PI times 2).
	constexpr real TAU = PI * 2;

	// Inverse of the PI constant.
	constexpr real INVPI = 1.0 / PI;

	// Square root of PI constant.
	constexpr real SQRTPI = 1.772454;

	// Euler constant.
	constexpr real E = 2.718281828459045235360287;

	// Logarithm base 2 of Euler constant.
	constexpr real LOG2E = 1.44269504088896338700465094;

	// Logarithm base 2 of 10.
	constexpr real LOG210 = 3.32192809488736218170856773213;

	// Logarithm base 10 of Euler constant.
	constexpr real LOG10E = 0.434294481903;

	// Natural logarithm of 2.
	constexpr real LN2 = 0.69314718056;

	// Natural logarithm of 10.
	constexpr real LN10 = 2.30258509299;

	// Multiply a number by this constant to change it to radians.
	constexpr real DEG2RAD = 0.017453292519943295474371680598;

	// Multiply a number by this constant to change it to degrees.
	constexpr real RAD2DEG = 57.2957795130823228646477218717;

	// Square root of 2.
	constexpr real SQRT2 = 1.41421356237;

	// Inverse square root of 2.
	constexpr real INVSQR2 = 0.707106781187;

	// Order of Taylor series approximations
	constexpr int TAYLOR_ORDER = UROBORO_TAYLOR_ORDER;

	// Default number of steps for integral approximation
	constexpr int INTEGRATION_STEPS = UROBORO_INTEGRATION_STEPS;

	// Relative precision for derivative approximation
	constexpr real DERIV_PREC = UROBORO_DERIV_PREC;

	// Biggest fractional part to ignore in powf computation
	constexpr real POW_APPROXIMATION_TOLERANCE = UROBORO_POWF_APPROX_TOL;

	// Approximation tolerance for root finding
	constexpr real ROOT_APPROX_TOL = UROBORO_ROOT_APPROX_TOL;

	// Approximation tolerance for bisection root finding
	constexpr real BISECTION_APPROX_TOL = UROBORO_BISECTION_APPROX_TOL;

	// Approximation tolerance for Newton-Raphson root finding
	constexpr real NEWTON_RAPHSON_TOL = UROBORO_NEWTON_RAPHSON_TOL;

	// Maximum number of iterations for bisection
	constexpr unsigned int MAX_BISECTION_ITER = UROBORO_MAX_BISECTION_ITER;

	// Maximum number of iterations for Newton-Raphson
	constexpr unsigned int MAX_NEWTON_ITER = UROBORO_MAX_NEWTON_ITER;

	// Maximum number of iterations for Steffensen
	constexpr unsigned int MAX_STEFFENSEN_ITER = UROBORO_MAX_STEFFENSEN_ITER;

	// Maximum number of iterations for Chebyshev
	constexpr unsigned int MAX_CHEBYSHEV_ITER = UROBORO_MAX_CHEBYSHEV_ITER;

}

#endif
