
///
/// @file autodiff_types.h Types and traits for automatic differentiation.
///

#ifndef THEORETICA_AUTODIFF_TYPES_H
#define THEORETICA_AUTODIFF_TYPES_H


namespace theoretica {


	// Type traits for automatic differentiation


	/// Type trait to check whether the given type is a multidual number
	template<typename Type>
	struct is_dual_type : std::false_type {};

	/// Type trait to check whether the given type is a multidual number
	template<>
	struct is_dual_type<dual> : std::true_type {};


	/// Type trait to check whether the given function takes
	/// a dual number as its first argument.
	template<typename Function>
	using is_dual_func =
	std::conditional_t <
		is_dual_type <
			typename _internal::return_type_or_void<Function, dual>::type
		>::value,
		std::true_type, std::false_type
	>;


	/// Enable a certain function overload if the given type
	/// is a function taking as first argument a dual number
	template<typename Function, typename T = bool>
	using enable_dual_func =
		typename std::enable_if<is_dual_func<Function>::value, T>::type;


	/// Type trait to check whether the given type is a multidual number
	template<typename Type>
	struct is_dual2_type : std::false_type {};

	/// Type trait to check whether the given type is a multidual number
	template<>
	struct is_dual2_type<dual2> : std::true_type {};


	/// Type trait to check whether the given function takes
	/// a dual2 number as its first argument.
	template<typename Function>
	using is_dual2_func =
	std::conditional_t <
		is_dual2_type <
			typename _internal::return_type_or_void<Function, dual2>::type
		>::value,
		std::true_type, std::false_type
	>;


	/// Enable a certain function overload if the given type
	/// is a function taking as first argument a dual2 number
	template<typename Function, typename T = bool>
	using enable_dual2_func =
		typename std::enable_if<is_dual2_func<Function>::value, T>::type;


	/// Type trait to check whether the given type is a multidual number
	template<typename Type>
	struct is_multidual_type : std::false_type {};

	/// Type trait to check whether the given type is a multidual number
	template<unsigned int N>
	struct is_multidual_type<multidual<N>> : std::true_type {};

	/// Type trait to check whether the given function takes a
	/// multidual vector as first argument.
	template<typename Function>
	using is_multidual_func =
		std::conditional_t<
			is_multidual_type<
				vector_element_or_void_t<
					std::tuple_element_t<
						0, typename extract_func_args<Function>::type
					>
				>
			>::value, std::true_type, std::false_type>;


	/// Enable a certain function overload if the given type
	/// is a function taking as first argument a multidual number
	template<typename Function, typename T = bool>
	using enable_multidual_func =
		typename std::enable_if<is_multidual_func<Function>::value, T>::type;


	namespace autodiff {


		// Alias types for multivariate automatic differentiation


		/// Real type for multivariate automatic differentiation
		/// (read "differential real").
		template<unsigned int N = 0>
		using dreal_t = multidual<N>;

		/// Vector type for multivariate automatic differentiation
		/// (read "differential vector").
		template<unsigned int N = 0>
		using dvec_t = vec<dreal_t<N>, N>;


		/// Real type for multivariate automatic differentiation
		/// with dynamically allocated vectors.
		using dreal = dreal_t<0>;

		/// Vector type for multivariate automatic differentiation
		/// with dynamically allocated vectors.
		using dvec = dvec_t<0>;


		/// Real type for multivariate automatic differentiation
		/// with two-dimensional statically allocated vectors.
		using dreal2 = dreal_t<2>;

		/// Vector type for multivariate automatic differentiation
		/// with two-dimensional statically allocated vectors.
		using dvec2 = dvec_t<2>;


		/// Real type for multivariate automatic differentiation
		/// with three-dimensional statically allocated vectors.
		using dreal3 = dreal_t<3>;

		/// Vector type for multivariate automatic differentiation
		/// with three-dimensional statically allocated vectors.
		using dvec3 = dvec_t<3>;


		/// Real type for multivariate automatic differentiation
		/// with four-dimensional statically allocated vectors.
		using dreal4 = dreal_t<4>;

		/// Vector type for multivariate automatic differentiation
		/// with four-dimensional statically allocated vectors.
		using dvec4 = dvec_t<4>;
	}
}

#endif
