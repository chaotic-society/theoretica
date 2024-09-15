
///
/// @file core/dataset.h Operations on datasets
///

#ifndef THEORETICA_DATASET_H
#define THEORETICA_DATASET_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "./core_traits.h"
#include "./constants.h"
#include "./error.h"


namespace theoretica {


	// Operations on datasets and generic ordered sets of numbers
	// The Vector type must have size() and operator[]() methods
	// (e.g. std::vector<real> and vec<real>)


	/// Compute the product of a set of values
	template<typename Vector, enable_vector<Vector> = true>
	inline auto product(const Vector& X) {

		if(!X.size()) {
			TH_MATH_ERROR("product", X.size(), INVALID_ARGUMENT);
			return vector_element_t<Vector>(nan());
		}

		auto res = X[0];
		for(unsigned int i = 1; i < X.size(); i++)
			res *= X[i];

		return res;
	}


	/// Sum the products of two sets of values
	template<typename Vector>
	inline auto product_sum(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size() || !X.size()) {
			TH_MATH_ERROR("product_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		auto res = X[0] * Y[0];
		for(unsigned int i = 1; i < X.size(); i++)
			res += X[i] * Y[i];

		return res;
	}


	/// Sum the products of the squares of two sets of data
	template<typename Vector>
	inline auto product_sum_squares(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size() || !X.size()) {
			TH_MATH_ERROR("product_sum_squares", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		auto res = (X[0] * X[0]) * (Y[0] * Y[0]);

		for(unsigned int i = 1; i < X.size(); i++)
			res += (X[i] * X[i]) * (Y[i] * Y[i]);

		return res;
	}


    /// Sum the products of three sets of values
    template<typename Vector>
	inline auto product_sum(const Vector& X, const Vector& Y, const Vector& Z) {

		if(X.size() != Y.size() || X.size() != Z.size() || !X.size()) {
			TH_MATH_ERROR("product_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		auto res = X[0] * Y[0] * Z[0];
		for(unsigned int i = 1; i < X.size(); i++)
			res += X[i] * Y[i] * Z[i];

		return res;
    }


    /// Sum the quotients of two sets of values
	template<typename Vector>
	inline auto quotient_sum(const Vector& X, const Vector& Y) {

		if(X.size() != Y.size()) {
			TH_MATH_ERROR("quotient_sum", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		if(abs(Y[0]) < MACH_EPSILON) {
			TH_MATH_ERROR("quotient_sum", Y[0], DIV_BY_ZERO);
			return nan();
		}

		auto res = X[0] / Y[0];
		for(unsigned int i = 1; i < X.size(); i++) {

			if(abs(Y[i]) < MACH_EPSILON) {
				TH_MATH_ERROR("quotient_sum", Y[i], DIV_BY_ZERO);
				return nan();
			}

			res += X[i] / Y[i];
		}

		return res;
	}


    /// Sum the squares of a set of values
    template<typename Vector>
	inline auto sum_squares(const Vector& X) {

		if(!X.size()) {
			TH_MATH_ERROR("sum_squares", X.size(), INVALID_ARGUMENT);
			return nan();
		}

		auto res = X[0] * X[0];
		for(unsigned int i = 1; i < X.size(); i++)
			res += X[i] * X[i];

		return res;
    }


    /// Compute the sum of a set of values using
	/// the compensated Neumaier-Kahan-Babushka summation
	/// algorithm to reduce round-off error.
	///
	/// @param X The vector of real values to sum
	template<typename Vector>
	inline real sum_compensated(const Vector& X) {

		// Total sum
		real sum = 0;

		// Correction term
		real corr = 0;

		for (unsigned int i = 0; i < X.size(); i++) {

			const real temp = sum + X[i];

			// Sort the two addends to preserve bits
			corr += (abs(sum) >= abs(X[i]))
			  ? ((sum - temp) + X[i])
			  : ((X[i] - temp) + sum);

			sum = temp;
		}

		// Apply correction
		return sum + corr;
	}


	/// Compute the sum of a set of values using
	/// pairwise summation to reduce round-off error.
	/// The function does not check for validity of begin and end indices.
	///
	/// @param X The vector of real values to sum
	/// @param begin The starting index of the sum, defaults to 0
	/// @param end One plus the last index of the sum, defaults to X.size()
	/// @param base_size The size of the base case, defaults to 128
	template<typename Vector>
	inline real sum_pairwise(
		const Vector& X, size_t begin = 0,
		size_t end = 0, size_t base_size = 128) {

		if(end == 0)
			end = X.size();

		real sum = 0;

		// Base case with given size (defaults to 128)
		if((end - begin) <= base_size) {

			for (size_t i = begin; i < end; ++i)
				sum += X[i];

		} else {

			// Recursive sum of two halves
			const size_t m = (end - begin) / 2;
			const size_t cutoff = begin + m;

			sum = sum_pairwise(X, begin, cutoff, base_size)
				+ sum_pairwise(X, cutoff, end, base_size);
		}

		return sum;
	}


	/// Compute the sum of a vector of real values
	/// using pairwise summation to reduce round-off error.
	///
	/// @param X The vector of values to sum
	template <
		typename Vector,
		std::enable_if_t<has_real_elements<Vector>::value> = true
	>
	inline auto sum(const Vector& X) {
		return sum_pairwise(X);
	}


	/// Compute the sum of a set of values.
	///
	/// @param X The vector of values to sum
	template<typename Vector>
	inline auto sum(const Vector& X) {

		// Use pairwise sum to reduce floating point errors.
		if TH_CONSTIF (has_real_elements<Vector>())
			return sum_pairwise(X);

		auto res = X[0];
		for(unsigned int i = 1; i < X.size(); i++)
			res += X[i];

		return res;
	}


	/// Apply a function to a set of values element-wise.
	/// @note Unlike functions in the vectorized namespace,
	/// this routine is not parallelized.
	///
	/// @param f The function to apply
	/// @param X The vector to modify the elements of
	/// @return A reference to the modified vector
	template<typename Vector, typename Function>
	inline Vector& apply(Function f, Vector& X) {

		for (unsigned int i = 0; i < X.size(); i++)
			X[i] = f(X[i]);

		return X;
	}


	/// Get a new vector obtained by applying the function element-wise.
	///
	/// @param f The function to apply
	/// @param src The input vector
	/// @param dest The output vector
	/// @return A reference to the modified output vector
	template<typename Vector1, typename Vector2 = Vector1, typename Function>
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


	/// Get a new vector obtained by applying the function element-wise.
	///
	/// @param f The function to apply
	/// @param X The input vector
	/// @return The modified output vector
	template<typename Vector1, typename Vector2 = Vector1, typename Function>
	inline Vector1 map(Function f, const Vector1& X) {

		Vector2 res;
		res.resize(X.size());

		for (unsigned int i = 0; i < X.size(); i++)
			res[i] = f(X[i]);

		return res;
	}


	/// Concatenate two datasets to form a single one
	template<typename Vector1, typename Vector2, typename Vector3 = Vector1>
	inline Vector3 concatenate(const Vector1& v1, const Vector2& v2) {

		Vector3 res;
		res.resize(v1.size() + v2.size());
		const unsigned int offset = v1.size();

		for (unsigned int i = 0; i < offset; ++i)
			res[i] = v1[i];

		for (unsigned int i = 0; i < v2.size(); ++i)
			res[i + offset] = v2[i];

		return res;
	}


	/// Finds the maximum value inside a dataset
	template<typename Vector>
	inline auto max(const Vector& X) {

		using Type = vector_element_t<Vector>;

		if(!X.size()) {
			TH_MATH_ERROR("max", X.size(), INVALID_ARGUMENT);
			return Type(nan());
		}

		auto curr = X[0];

		for (unsigned int i = 1; i < X.size(); ++i)
			if(X[i] > curr)
				curr = X[i];

		return curr;
	}


	/// Finds the minimum value inside a dataset
	template<typename Vector>
	inline auto min(const Vector& X) {

		using Type = vector_element_t<Vector>;

		if(!X.size()) {
			TH_MATH_ERROR("min", X.size(), INVALID_ARGUMENT);
			return Type(nan());
		}

		auto curr = X[0];

		for (unsigned int i = 1; i < X.size(); ++i)
			if(X[i] < curr)
				curr = X[i];

		return curr;
	}


	// Different types of means


	/// Compute the arithmetic mean of a set of values
	template<typename Dataset>
	inline real arithmetic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("arithmetic_mean", data.size(), DIV_BY_ZERO);
			return nan();
		}

		// Sum of x_i / N
		return sum(data) / (real) data.size();
	}


	/// Compute the harmonic mean of a set of values
	template<typename Dataset>
	inline real harmonic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("harmonic_mean", data.size(), DIV_BY_ZERO);
			return nan();
		}

		real res = 0;

		for (unsigned int i = 0; i < data.size(); ++i) {

			if(data[i] == 0) {
				TH_MATH_ERROR("harmonic_mean", data[i], DIV_BY_ZERO);
				return nan();
			}

			res += 1.0 / data[i];
		}

		return static_cast<real>(data.size()) / res;
	}


	/// Compute the geometric mean of a set of values
	/// as \f$\sqrt[n]{\Pi_i x_i}\f$
	template<typename Dataset>
	inline real geometric_mean(const Dataset& data) {
		return root(product(data), data.size());
	}


	/// Compute the weighted mean of a set of values
	/// <data> and <weights> must have the same size
	template<typename Dataset1, typename Dataset2>
	inline real weighted_mean(const Dataset1& data, const Dataset2& weights) {

		// Sum of x_i * w_i / Sum of w_i
		return product_sum(data, weights) / sum(weights);
	}


	/// Compute the quadratic mean (Root Mean Square) of a set of values
	/// \f$m_q = \sqrt{x1^2 + x2^2 + ...}\f$
	template<typename Dataset>
	inline real quadratic_mean(const Dataset& data) {

		if(!data.size()) {
			TH_MATH_ERROR("quadratic_mean", data.size(), INVALID_ARGUMENT);
			return nan();
		}

		return sqrt(sum_squares(data) / data.size());
	}
}

#endif
