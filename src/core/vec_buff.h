
///
/// @file core/vec_buff.h Operations on datasets
///

#ifndef THEORETICA_VEC_BUFF_H
#define THEORETICA_VEC_BUFF_H

#include <vector>
#include "./constants.h"
#include <functional>

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif


namespace theoretica {


	/// A dynamically allocated variable-size container.
	/// Defined by default as an alias for `std::vector<real>`
	using vec_buff = std::vector<real>;


	// Operations on datasets
	// The Vector type must have size() and operator[]() functions
	// (e.g. vec_buff and vec<real, N>)


	/// Compute the product of a set of values
	template<typename Vector>
	inline real product(const Vector& X) {

		if(!X.size()) {
			TH_MATH_ERROR("product", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 1;
		for(unsigned int i = 0; i < X.size(); i++)
			res *= X[i];

		return res;
	}


	/// Sum the products of two sets of values
	template<typename Vector>
	inline real product_sum(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("product_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++)
			res += X[i] * Y[i];

		return res;
	}


	/// Sum the products of the squares of two sets of data
	template<typename Vector>
	inline real product_sum_squares(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("product_sum_squares", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++)
			res += square(X[i]) * square(Y[i]);

		return res;
	}


    /// Sum the products of three sets of values
    template<typename Vector>
	inline real product_sum(const Vector& X, const Vector& Y, const Vector& Z) {

		if(X.size() != Y.size() || X.size() != Z.size()) {
			TH_MATH_ERROR("product_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++)
			res += X[i] * Y[i] * Z[i];

		return res;
    }


    /// Sum the quotients of two sets of values
	template<typename Vector>
	inline real quotient_sum(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("quotient_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++) {

			if(Y[i] == 0) {
				TH_MATH_ERROR("quotient_sum", Y[i], DIV_BY_ZERO);
				return nan();
			}

			res += X[i] / Y[i];
		}

		return res;
	}


    /// Sum the squares of a set of values
    template<typename Vector>
	inline real sum_squares(const Vector& X) {

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++)
			res += X[i] * X[i];

		return res;
    }


    /// Sum together a set of values
	template<typename Vector>
	inline real sum(const Vector& X) {

		real res = 0;
		for(unsigned int i = 0; i < X.size(); i++)
			res += X[i];

		return res;
	}


	/// Apply a function to a set of values element-wise
	/// (vectorized evaluation)
	template<typename Function, typename Vector>
	inline Vector& apply(Function f, Vector& X) {

		for (unsigned int i = 0; i < X.size(); i++)
			X[i] = f(X[i]);

		return X;
	}


	/// Get a new vector obtained by applying
	/// the function element-wise
	/// (vectorized evaluation)
	template<typename Function, typename Vector1, typename Vector2 = Vector1>
	inline Vector2& map(Function f, const Vector1& src, Vector2& dest) {

		if(src.size() != dest.size()) {
			TH_MATH_ERROR("th::map", dest.size(), INVALID_ARGUMENT);
			dest = Vector2(nan());
			return dest;
		}

		for (unsigned int i = 0; i < src.size(); i++)
			dest[i] = f(src[i]);

		return dest;
	}


	/// Get a new vector obtained by applying
	/// the function element-wise
	/// (vectorized evaluation)
	template<typename Function, typename Vector1, typename Vector2 = Vector1>
	inline Vector1 map(Function f, const Vector1& X) {

		Vector2 res;
		res.resize(X.size());

		for (unsigned int i = 0; i < X.size(); i++)
			res[i] = f(X[i]);

		return res;
	}


	/// Finds the maximum value inside a dataset
	template<typename Vector>
	inline real max(const Vector& X) {

		if(!X.size()) {
			TH_MATH_ERROR("max", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real curr = X[0];

		for (unsigned int i = 1; i < X.size(); ++i)
			if(X[i] > curr)
				curr = X[i];

		return curr;
	}


	/// Finds the minimum value inside a dataset
	template<typename Vector>
	inline real min(const Vector& X) {

		if(!X.size()) {
			TH_MATH_ERROR("min", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		real curr = X[0];

		for (unsigned int i = 1; i < X.size(); ++i)
			if(X[i] < curr)
				curr = X[i];

		return curr;
	}


#ifndef THEORETICA_NO_PRINT

	/// Stream the buffer in string representation to an output stream (std::ostream)
	inline std::ostream& operator<<(std::ostream& out, const std::vector<real>& obj) {

		for (real x : obj) {
			out << x << std::endl;
		}

		return out;
	}


	/// Stream a vector in string representation to an output stream (std::ostream)
	// template<typename T>
	// inline std::ostream& operator<<(std::ostream& out, const std::vector<T>& obj) {

	// 	for (T x : obj) {
	// 		out << x << std::endl;
	// 	}

	// 	return out;
	// }

#endif

}

#endif
