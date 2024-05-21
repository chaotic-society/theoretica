
///
/// @file core/function.h Mathematical function pointers
///

#ifndef THEORETICA_FUNCTION_H
#define THEORETICA_FUNCTION_H

#include "./constants.h"
#include "../complex/complex.h"
#include "../algebra/vec.h"

#include <functional>


namespace theoretica {


	/// Function pointer to a real function of real variable
	using real_function = std::function<real(real)>;


	/// Function pointer to a complex function of complex variable
	using complex_function = std::function<complex<>(complex<>)>;


	/// Function pointer to a probability distribution function
	/// where the first argument is the variable and the second
	/// argument is a vector of the parameters of the distribution.
	using stat_function = std::function<real(real, const vec<real>&)>;

}



#endif
