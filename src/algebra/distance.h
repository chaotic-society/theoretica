
///
/// @file distance.h Distances and norms
///

#ifndef THEORETICA_DISTANCE_H
#define THEORETICA_DISTANCE_H

#include "../core/constants.h"
#include "../core/core_traits.h"
#include "../core/real_analysis.h"
#include "../complex/complex.h"
#include "../complex/complex_analysis.h"
#include "./vec.h"


namespace theoretica {

	/// Norms

	/// Compute the Lp norm of a vector.
	/// \f$L_p(\vec v) = (\Sigma_i \ |v_i|^p)^{1/p}\f$
	template<typename Vector>
	inline real lp_norm(const Vector& v, unsigned int p) {

		real sum = 0;

		for (unsigned int i = 0; i < v.size(); ++i)
			sum += pow(abs(v[i]), p);

		return root(sum, p);
	}


	/// Compute the L1 norm of a vector.
	/// \f$L_1(\vec v) = \Sigma_i \ |v_i|\f$
	template<typename Vector>
	inline real l1_norm(const Vector& v) {

		real sum = 0;

		for (unsigned int i = 0; i < v.size(); ++i)
			sum += abs(v[i]);

		return sum;
	}

	/// Compute the L2 norm of a vector.
	/// \f$L_2(\vec v) = \sqrt{\Sigma_i \ v_i^2}\f$
	template<typename Vector>
	inline real l2_norm(const Vector& v) {

		using Type = indexable_element_t<Vector>;
		Type sum = 0;

		if TH_CONSTIF (!is_complex_type<Type>())
			for (unsigned int i = 0; i < v.size(); ++i)
				sum += v[i] * conjugate(v[i]);
		else
			for (unsigned int i = 0; i < v.size(); ++i)
				sum += square(v[i]);

		return sqrt(sum);
	}


	/// Compute the Linf norm of a vector.
	/// \f$L_{inf}(\vec v) = max(|v_i|)\f$
	template<typename Vector>
	inline real linf_norm(const Vector& v) {

		real res = abs(v[0]);

		for (unsigned int i = 1; i < v.size(); ++i)
			res = max(res, abs(v[i]));

		return res;
	}


	/// Distances

	/// Compute the Euclidean distance between two vectors
	template<typename Vector>
	inline real euclidean_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("euclidean_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		return l2_norm(v1 - v2);
	}


	/// Compute the Euclidean distance between two vectors
	template<unsigned int N>
	inline real distance(const vec<real, N>& v1, const vec<real, N>& v2) {
		return euclidean_distance(v1, v2);
	}


	/// Compute the Euclidian distance between two values
	inline real euclidean_distance(real a, real b) {
		return abs(a - b);
	}


	/// Compute the Euclidian distance between two values
	inline real distance(real a, real b) {
		return euclidean_distance(a, b);
	}


	/// Compute the distance between two complex numbers
	template<typename T>
	inline complex<T> distance(complex<T> z1, complex<T> z2) {
		return (z1 - z2).norm();
	}


	/// Compute the Minkowski distance between two vectors
	template<typename Vector>
	inline real minkowski_distance(const Vector& v1, const Vector& v2, unsigned int p) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("minkowski_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		return lp_norm(v1 - v2, p);
	}


	/// Compute the Minkowski distance between two values
	inline real minkowski_distance(real a, real b, unsigned int p) {
		return root(pow(abs(b - a), p), p);
	}


	/// Compute the Hermitian distance between two vectors
	template<typename Vector, typename T = real>
	inline auto hermitian_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("hermitian_distance", v1.size(), INVALID_ARGUMENT);
			return  complex<T>(nan());
		}

		complex<T> sum = 0;

		for (size_t i = 0; i < v1.size(); ++i) {
			const complex<T> diff = v1[i] - conjugate(v2[i]);
			sum += diff * diff.conjugate();
		}

		return sqrt(sum);
	}


	/// Compute the Hermitian distance between two vectors
	template<unsigned int N, typename T>
	inline complex<T> distance(
		const vec<complex<T>, N>& v1, const vec<complex<T>, N>& v2) {

		return hermitian_distance<vec<complex<T>, N>, T>(v1, v2);
	}


	/// Compute the Manhattan distance between two vectors
	template<typename Vector>
	inline real manhattan_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("manhattan_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		return l1_norm(v1 - v2);
	}


	/// Compute the Chebyshev distance between two vectors
	template<typename Vector>
	inline real chebyshev_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("chebyshev_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		return linf_norm(v1 - v2);
	}


	/// Compute the discrete distance between two vectors
	template<typename Vector>
	inline real discrete_distance(
		const Vector& v1, const Vector& v2, real tolerance = MACH_EPSILON) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("discrete_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		bool diff = false;

		for (size_t i = 0; i < v1.size(); ++i) {

			// The vectors differ if a pair of elements differs
			// more than the given tolerance
			if(abs(v1[i] - v2[i]) > tolerance) {
				diff = true;
				break;
			}
		}

		return diff ? 1 : 0;
	}


	/// Compute the Canberra distance between two vectors
	template<typename Vector>
	inline real canberra_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("canberra_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		real sum = 0;

		for (size_t i = 0; i < v1.size(); ++i)
			sum += abs(v1[i] - v2[i]) / (abs(v1[i]) + abs(v2[i]));

		return sum;
	}


	/// Compute the cosine distance between two vectors
	template<typename Vector>
	inline real cosine_distance(const Vector& v1, const Vector& v2) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("cosine_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		real product = 0;
		real sum_sqr_x = 0;
		real sum_sqr_y = 0;

		for (size_t i = 0; i < v1.size(); ++i) {
			product += v1[i] * v2[i];
			sum_sqr_x += square(v1[i]);
			sum_sqr_y += square(v2[i]);
		}

		return product / sqrt(sum_sqr_x * sum_sqr_y);
	}


	/// Compute the Hamming distance between two vectors
	template<typename Vector>
	inline real hamming_distance(
		const Vector& v1, const Vector& v2, real tolerance = MACH_EPSILON) {

		if(v1.size() != v2.size()) {
			TH_MATH_ERROR("hamming_distance", v1.size(), INVALID_ARGUMENT);
			return nan();
		}

		unsigned int count = 0;

		// Count how many elements differ to some tolerance
		for (unsigned int i = 0; i < v1.size(); ++i)
			if(abs(v1[i] - v2[i]) > tolerance)
				count++;

		return count;
	}

}

#endif
