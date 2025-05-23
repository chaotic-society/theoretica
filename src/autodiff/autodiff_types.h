
///
/// @file autodiff_types.h Types and traits for automatic differentiation.
///

#ifndef THEORETICA_AUTODIFF_TYPES_H
#define THEORETICA_AUTODIFF_TYPES_H


namespace theoretica {

	namespace autodiff {


		// Type traits for automatic differentiation


		// Type trait to check whether the given type is a multidual number
		template<typename Type>
		struct is_dual_type : std::false_type {};

		
		template<>
		struct is_dual_type<dual> : std::true_type {};


		// Type trait to check whether the given function takes
		// a dual number as its first argument.
		template<typename Function>
		using is_dual_func =
		std::conditional_t <
			is_dual_type <
				typename _internal::return_type_or_void<Function, dual>::type
			>::value,
			std::true_type, std::false_type
		>;


		// Enable a certain function overload if the given type
		// is a function taking as first argument a dual number
		template<typename Function, typename T = bool>
		using enable_dual_func =
			typename std::enable_if<is_dual_func<Function>::value, T>::type;


		// Type trait to check whether the given type is a multidual number
		template<typename Type>
		struct is_dual2_type : std::false_type {};

		
		template<>
		struct is_dual2_type<dual2> : std::true_type {};


		// Type trait to check whether the given function takes
		// a dual2 number as its first argument.
		template<typename Function>
		using is_dual2_func =
		std::conditional_t <
			is_dual2_type <
				typename _internal::return_type_or_void<Function, dual2>::type
			>::value,
			std::true_type, std::false_type
		>;


		// Enable a certain function overload if the given type
		// is a function taking as first argument a dual2 number
		template<typename Function, typename T = bool>
		using enable_dual2_func =
			typename std::enable_if<is_dual2_func<Function>::value, T>::type;


		// Type trait to check whether the given type is a multidual number
		template<typename Type>
		struct is_multidual_type : std::false_type {};

		
		template<unsigned int N>
		struct is_multidual_type<multidual<N>> : std::true_type {};

		// Enable a certain function overload if the given type
		// is an instantiation of the multidual template class
		template<typename Type, typename T = bool>
		using enable_multidual = 
			typename std::enable_if<is_multidual_type<Type>::value, T>::type;


		// Enable a certain function overload if the given type
		// is a Callable object corresponding to a multidual function
		// representing a scalar field, that is, a function taking
		// a dvec_t and returning a dreal_t.
		template<typename Function, typename T = bool>
		using enable_scalar_field = typename std::enable_if <
			is_multidual_type<return_type_t<Function>>::value, T
		>::type;


		// Enable a certain function overload if the given type
		// is a Callable object corresponding to a multidual function
		// representing a vector field, that is, a function taking
		// a dvec_t and returning a dvec_t.
		template<typename Function, typename T = bool>
		using enable_vector_field = typename std::enable_if <
			is_multidual_type <
				vector_element_t<return_type_t<Function>>
			>::value, T
		>::type;


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
