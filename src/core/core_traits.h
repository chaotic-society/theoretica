
///
/// @file core_traits.h Fundamental type traits.
///

#ifndef THEORETICA_CORE_TRAITS_H
#define THEORETICA_CORE_TRAITS_H

#include <type_traits>
#include <tuple>
#include "constants.h"


namespace theoretica {


	// Implement a simple void_t trait in C++14.
	/// @internal
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


	/// Check whether a structure is orderable,
	/// by checking that it has a comparison operator<().
	template<typename Structure, typename = _internal::void_t<>>
	struct is_orderable : std::false_type{};

	// The @internal directive is used to avoid documenting this overload
	/// @internal
	template<typename Structure>
	struct is_orderable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>() < std::declval<Structure>())>
	> : std::true_type{};


	/// Check whether a structure is indexable by a single integer index,
	/// by checking that it has the operator[](0).
	template<typename Structure, typename = _internal::void_t<>>
	struct is_indexable : std::false_type{};

	/// @internal
	template<typename Structure>
	struct is_indexable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>()[0])>
	> : std::true_type{};


	/// Check whether a structure is iterable,
	/// by checking that it has a method begin().
	template<typename Structure, typename = _internal::void_t<>>
	struct is_iterable : std::false_type{};

	/// @internal
	template<typename Structure>
	struct is_iterable
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>().begin())>
	> : std::true_type{};


	/// Check whether a structure is considerable a vector,
	/// by checking that it has an operator[] and a size() method.
	template<typename Structure, typename = _internal::void_t<>>
	struct is_vector : std::false_type{};

	template<typename Structure>
	struct is_vector
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>()[0]),
		decltype(std::declval<Structure>().size())>
	> : std::true_type{};


	/// Check whether a structure is considerable a matrix,
	/// by checking that it has an operator(), a rows() method
	/// and a cols() method.
	template<typename Structure, typename = _internal::void_t<>>
	struct is_matrix : std::false_type{};

	/// @internal
	template<typename Structure>
	struct is_matrix
	<Structure, _internal::void_t<
		decltype(std::declval<Structure>()(0, 0)),
		decltype(std::declval<Structure>().rows()),
		decltype(std::declval<Structure>().cols())>
	> : std::true_type{};


	/// @internal
	namespace _internal {

		/// Helper structure for vector_element_t
		template<typename Structure, typename = void>
		struct vector_element_or_void {
			using type = void;
		};

		/// Helper structure for vector_element_t
		template<typename Structure>
		struct vector_element_or_void
			<Structure, _internal::void_t<decltype(std::declval<Structure&>()[0])>> {
			using type = std::remove_reference_t<decltype(std::declval<Structure>()[0])>;
		};
	}

	/// Extract the type of a vector (or any indexable container) from its operator[],
	/// returning void if the type has no operator[].
	template<typename Structure>
	using vector_element_or_void_t =
		typename _internal::vector_element_or_void<Structure>::type;


	/// Extract the type of a vector (or any indexable container) from its operator[].
	template<typename Structure>
	using vector_element_t =
		std::remove_reference_t<decltype(std::declval<Structure>()[0])>;


	/// Extract the type of a matrix (or any doubly indexable container) from its operator().
	template<typename Structure>
	using matrix_element_t =
		std::remove_reference_t<decltype(std::declval<Structure>()(0, 0))>;


	/// Type trait to check whether an indexable container
	/// has elements of the given type.
	template<typename Structure, typename Type>
	struct has_type_elements
	: std::is_same<vector_element_t<Structure>, Type> {};


	/// Type trait to check whether an indexable container
	/// has real elements.
	template<typename Structure>
	using has_real_elements = is_real_type<vector_element_t<Structure>>;


	/// Enable a function overload if the template typename
	/// is considerable a matrix. The std::enable_if structure
	/// is used, with type T which defaults to bool.
	template<typename Structure, typename T = bool>
	using enable_matrix = std::enable_if_t<is_matrix<Structure>::value, T>;


	/// Enable a function overload if the template typename
	/// is considerable a vector. The std::enable_if structure
	/// is used, with type T which defaults to bool.
	template<typename Structure, typename T = bool>
	using enable_vector = std::enable_if_t<is_vector<Structure>::value, T>;


	/// Extract the type of the arguments of a function.
	template<typename Function>
	struct extract_func_args;

	/// @internal
	template<typename Function, typename... Args>
	struct extract_func_args<Function(Args...)> {
		using type = std::tuple<Args...>;
	};


	namespace _internal {

		// Helper structure to extract the first type of a variadic template
		template<typename Arg, typename ...Other>
		struct get_first {
			using type = Arg;
		};

		// Helper structure for is_real_func and related traits
		template<typename Function, typename T, typename = void>
		struct return_type_or_void {
			using type = void;
		};

		// Helper structure for is_real_func and related traits
		template<typename Function, typename T>
		struct return_type_or_void
			<Function, T, _internal::void_t<decltype(std::declval<Function>()(T(0.0)))>> {
			using type = decltype(std::declval<Function>()(T(0.0)));
		};


		// Helper structure for extracting information from Callables
		template<typename T>
		struct func_helper : public func_helper<decltype(&T::operator())> {};

		template<typename ReturnType, typename ...Args>
		struct func_helper<ReturnType(Args...)> {

			using return_type = ReturnType;
			using args_type = std::tuple<Args...>;
			using first_arg_type = typename get_first<Args...>::type;
		};

		template<typename ReturnType, typename ...Args>
		struct func_helper<ReturnType(*)(Args...)> {

			using return_type = ReturnType;
			using args_type = std::tuple<Args...>;
			using first_arg_type = typename get_first<Args...>::type;
		};

		template<typename Class, typename ReturnType, typename ...Args>
		struct func_helper<ReturnType(Class::*)(Args...) const> {

			using return_type = ReturnType;
			using args_type = std::tuple<Args...>;
			using first_arg_type = typename get_first<Args...>::type;
		};
	}


	/// Type trait to check whether the given function takes
	/// a real number as its first argument.
	template<typename Function>
	using is_real_func =
	std::conditional_t<
		is_real_type <
			typename _internal::return_type_or_void<Function, real>::type
		>::value,
		std::true_type, std::false_type
	>;


	/// Enable a certain function overload if the given type
	/// is a function taking as first argument a real number
	template<typename Function, typename T = bool>
	using enable_real_func =
		typename std::enable_if_t<is_real_func<Function>::value, T>;


	/// Extract the return type of a Callable object, such as a function
	/// pointer or lambda function.
	template<typename Function>
	using return_type_t = typename _internal::func_helper<Function>::return_type;

}

#endif
