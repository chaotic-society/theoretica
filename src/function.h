
///
/// @file function.h Mathematical function pointers
///

#ifndef UROBORO_FUNCTION_H
#define UROBORO_FUNCTION_H

#include "./constants.h"
#include "./complex.h"


namespace uroboro {

	/// Real function
	using real_function = real(*)(real);

	/// Complex function
	using complex_function = complex(*)(complex);

	/// Multi-variable function
	using multi_function = real(*)(const vec_buff&);

	/// Statistical function
	using stat_function = real(*)(real, const vec_buff&);

}



#endif
