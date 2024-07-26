
///
/// @file core_traits.h Fundamental type traits.
///

#ifndef THEORETICA_CORE_TRAITS_H
#define THEORETICA_CORE_TRAITS_H

#include <type_traits>
#include "constants.h"


namespace theoretica {


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
	using indexable_element_t =
		std::remove_reference_t<decltype(std::declval<Structure>()[0])>;


	/// Extract the type of a matrix (or any doubly indexable container)
	/// from its operator().
	template<typename Structure>
	using matrix_element_t =
		std::remove_reference_t<decltype(std::declval<Structure>()(0, 0))>;


	/// Type trait to check whether an indexable container
	/// has elements of the given type.
	template<typename Structure, typename Type>
	struct has_type_elements
	: std::is_same<indexable_element_t<Structure>, Type> {};


	/// Type trait to check whether an indexable container
	/// has real elements.
	template<typename Structure>
	using has_real_elements = is_real_type<indexable_element_t<Structure>>;


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

}

#endif
