
///
/// @file parallel.h Parallelized methods to evaluate a function over a vector element-wise.
///

#ifndef THEORETICA_PARALLEL_H
#define THEORETICA_PARALLEL_H

#include "./vec.h"
#include "../core/dataset.h"
#include "../core/real_analysis.h"
#include "../complex/complex_analysis.h"
#include "../autodiff/dual_functions.h"
#include "../autodiff/dual2_functions.h"
#include "../autodiff/multidual_functions.h"


namespace theoretica {


	/// @namespace theoretica::parallel Parallelized element-wise evaluation of functions.
	namespace parallel {


		// Unary operations


		/// Parallel element-wise evaluation of a function,
		/// using OpenMP to speed up execution over a vector.
		///
		/// @param f The function to evaluate
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Function, typename Vector>
		inline Vector apply_function(Function f, const Vector& v) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = f(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the square function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector square(const Vector& v) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::square(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the cube function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector cube(const Vector& v) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::cube(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the abs function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector abs(const Vector& v) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::abs(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the pow function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector pow(const Vector& v, int n) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::pow(v[i], n);

			return res;
		}


		/// Parallel element-wise evaluation of the powf function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector powf(const Vector& v, real x) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::powf(v[i], x);

			return res;
		}


		/// Parallel element-wise evaluation of the sqrt function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector sqrt(const Vector& v) {

			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::sqrt(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the cbrt function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector cbrt(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::cbrt(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the exp function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector exp(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::exp(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the ln function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector ln(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::ln(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the log2 function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector log2(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::log2(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the log10 function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector log10(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::log10(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the sin function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector sin(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::sin(v[i]);

			return res;
		}


		/// Parallel element-wise evaluation of the cos function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector cos(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::cos(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the tan function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector tan(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::tan(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the cot function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector cot(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::cot(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the asin function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector asin(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::asin(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the acos function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector acos(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::acos(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the atan function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector atan(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::atan(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the sinh function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector sinh(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::sinh(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the cosh function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector cosh(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::cosh(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the tanh function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector tanh(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::tanh(v[i]);

			return res;
		}


		/// Vectorized (element-wise evaluation) of the coth function.
		///
		/// @param v The vector of inputs
		/// @return The transformed vector
		template<typename Vector>
		inline Vector coth(const Vector& v) {
			
			Vector res;
			res.resize(v.size());

			#pragma omp parallel for
			for (unsigned int i = 0; i < v.size(); i++)
				res[i] = theoretica::coth(v[i]);

			return res;
		}
	}
}

#endif
