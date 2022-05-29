
///
/// @file core/function.h Mathematical function pointers
///

#ifndef THEORETICA_FUNCTION_H
#define THEORETICA_FUNCTION_H

#include "./constants.h"
#include "../complex/complex.h"
#include "./vec_buff.h"


namespace theoretica {

	/// Function pointer to a real function
	using real_function = real(*)(real);

	/// Function pointer to a complex function
	using complex_function = complex(*)(complex);

	/// Function pointer to a statistical function
	using stat_function = real(*)(real, const vec_buff&);

}



#endif
