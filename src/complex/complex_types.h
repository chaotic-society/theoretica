
///
/// @file complex_types.h Complex data types definitions
///

#ifndef THEORETICA_COMPLEX_TYPES_H
#define THEORETICA_COMPLEX_TYPES_H

#include <type_traits>
#include "../core/core_traits.h"
#include "../core/constants.h"
#include "./complex.h"


namespace theoretica {


	/// A bi-complex number
	template<typename Type = real>
	using bicomplex = complex<complex<Type>>;


	/// Type trait to check whether the given type is a specialization
	/// of the complex number class or not, using the
	/// static boolean element is_complex_type<T>::value
	template<typename T>
	struct is_complex_type : std::false_type {};

	/// Type trait to check whether the given type is a specialization
	/// of the complex number class or not, using the
	/// static boolean element is_complex_type<T>::value
	template<typename T>
	struct is_complex_type<complex<T>> : std::true_type {};


	/// Type trait to check whether a container has complex elements.
	template<typename Structure>
	using has_complex_elements = is_complex_type<indexable_element_t<Structure>>;

}


#endif
