
///
/// @file vec_functions.h Parallelized methods to evaluate a function over a vector element-wise.
///

#ifndef THEORETICA_PARALLEL_VEC_FUNCTIONS_H
#define THEORETICA_PARALLEL_VEC_FUNCTIONS_H

#include "core/real_analysis.h"
#include "complex/complex_analysis.h"
#include "autodiff/dual_functions.h"
#include "autodiff/dual2_functions.h"
#include "autodiff/multidual_functions.h"


namespace theoretica {


/// @namespace theoretica::parallel Parallelized element-wise evaluation of functions.
namespace parallel {


	/// Parallel element-wise evaluation of a function,
	/// using OpenMP to speed up execution over a vector.
	///
	/// @param f The function to evaluate
	/// @param v The vector of inputs
	/// @return The transformed vector
	template<typename Function, typename Vector>
	inline Vector map(const Vector& v, Function f) {

		Vector res;
		res.resize(v.size());

		#pragma omp parallel for
		for (unsigned int i = 0; i < v.size(); i++)
			res[i] = f(v[i]);

		return res;
	}


	/// Overwrite a vector by applying a function to the elements
	/// of another vector. The two vectors must have the same size.
	///
	/// @param f The function to apply
	/// @param src The input vector
	/// @param dest The output vector
	/// @return A reference to the modified output vector
	template<typename Vector1, typename Vector2 = Vector1, typename Function>
	inline Vector2& map(const Vector1& src, Vector2& dest, Function f) {

		if(src.size() != dest.size()) {
			TH_MATH_ERROR("th::map", dest.size(), MathError::InvalidArgument);
			return algebra::vec_error(dest);
		}

		#pragma omp parallel for
		for (unsigned int i = 0; i < src.size(); i++)
			dest[i] = f(src[i]);

		return dest;
	}


	/// Transform a vector by applying a given function
	/// to each of its elements, modifying the vector itself
	/// and using OpenMP to speed up execution.
	///
	/// @param f The function to evaluate
	/// @param v The vector of inputs
	/// @return A reference to the input vector which has been modified
	template<typename Function, typename Vector>
	inline Vector& transform(Vector& v, Function f) {

		#pragma omp parallel for
		for (unsigned int i = 0; i < v.size(); i++)
			v[i] = f(v[i]);

		return v;
	}

}}

#endif
