
///
/// @file vec_functions.h Vectorized methods to evaluate a function over a vector element-wise.
///

#ifndef THEORETICA_VEC_FUNCTIONS_H
#define THEORETICA_VEC_FUNCTIONS_H

#include "../core/core_traits.h"
#include "../core/real_analysis.h"
#include "../complex/complex_analysis.h"
#include "../autodiff/dual_functions.h"
#include "../autodiff/dual2_functions.h"
#include "../autodiff/multidual_functions.h"

#ifdef _OPENMP
#include "../parallel/algebra/vec_functions.h"
#endif


namespace theoretica {


    // Unary operations


    /// Parallel element-wise evaluation of a function,
    /// using OpenMP to speed up execution over a vector.
    ///
    /// @param f The function to evaluate
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Function,
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector map(const Vector& v, Function f) {

        // Parallel dispatch
        #ifdef _OPENMP
        if (v.size() >= 1000)
            return parallel::map(v, f);
        #endif

        Vector res;
        res.resize(v.size());

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

        // Parallel dispatch
        #ifdef _OPENMP
        if (src.size() >= 1000)
            return parallel::map(src, dest, f);
        #endif

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
    template <
        typename Function,
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector& transform(Vector& v, Function f) {

        // Parallel dispatch
        #ifdef _OPENMP
        if (v.size() >= 1000)
            return parallel::transform(v, f);
        #endif

        for (unsigned int i = 0; i < v.size(); i++)
            v[i] = f(v[i]);

        return v;
    }


    /// Parallel element-wise evaluation of the square function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector square(const Vector& v) {

        return map(
            v, [](const auto& x) { return square(x); }
        );
    }


    /// Parallel element-wise evaluation of the cube function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector cube(const Vector& v) {

        return map(
            v, [](const auto& x) { return cube(x); }
        );
    }


    /// Parallel element-wise evaluation of the abs function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector abs(const Vector& v) {

        return map(
            v, [](const auto& x) { return abs(x); }
        );
    }


    /// Parallel element-wise evaluation of the pow function.
    ///
    /// @param v The vector of inputs
    /// @param n The integer exponent
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector pow(const Vector& v, int n) {

        return map(
            v, [n](const auto& x) { return pow(x, n); }
        );
    }


    /// Parallel element-wise evaluation of the powf function.
    ///
    /// @param v The vector of inputs
    /// @param x The real exponent
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector powf(const Vector& v, real x) {

        return map(
            v, [x](const auto& y) { return theoretica::powf(y, x); }
        );
    }


    /// Parallel element-wise evaluation of the sqrt function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector sqrt(const Vector& v) {

        return map(
            v, [](const auto& x) { return sqrt(x); }
        );
    }


    /// Parallel element-wise evaluation of the cbrt function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector cbrt(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return cbrt(x); }
        );
    }


    /// Parallel element-wise evaluation of the exp function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector exp(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return exp(x); }
        );
    }


    /// Parallel element-wise evaluation of the ln function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector ln(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return ln(x); }
        );
    }


    /// Parallel element-wise evaluation of the log2 function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector log2(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return log2(x); }
        );
    }


    /// Parallel element-wise evaluation of the log10 function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector log10(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return log10(x); }
        );
    }


    /// Parallel element-wise evaluation of the sin function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector sin(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return sin(x); }
        );
    }


    /// Parallel element-wise evaluation of the cos function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector cos(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return cos(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the tan function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector tan(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return tan(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the cot function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector cot(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return cot(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the asin function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector asin(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return asin(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the acos function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector acos(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return acos(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the atan function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector atan(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return atan(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the sinh function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector sinh(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return sinh(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the cosh function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector cosh(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return cosh(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the tanh function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector tanh(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return tanh(x); }
        );
    }


    /// Vectorized (element-wise evaluation) of the coth function.
    ///
    /// @param v The vector of inputs
    /// @return The transformed vector
    template <
        typename Vector,
        enable_vector<Vector> = true
    >
    inline Vector coth(const Vector& v) {
        
        return map(
            v, [](const auto& x) { return coth(x); }
        );
    }
}

#endif
