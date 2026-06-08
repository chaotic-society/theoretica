
///
/// @file quasirandom.h Quasi-random sequences
///

#ifndef THEORETICA_QUASIRANDOM_H
#define THEORETICA_QUASIRANDOM_H

#include "../core/real_analysis.h"
#include "../algebra/algebra_types.h"
#include "../optimization/roots.h"


namespace theoretica {
namespace random {


	/// Weyl quasi-random sequence in 1 dimension.
	///
	/// @param n The index of the element in the sequence
	/// @param alpha The base of the sequence, defaults to 
	/// the inverse of the Golden Section.
	/// @param shift The amount to shift the sequence (defaults to 0)
	/// 
	/// The Weyl sequence is defined as \f$x_n = \{n \alpha\}\f$,
	/// where \f$\{ \}\f$ is the fractional part.
	/// @note The alpha argument should be an irrational number.
	inline real weyl(unsigned int n, real alpha = INVPHI, real shift = 0.0) {
		return fract(n * (alpha + shift));
	}


	/// @class WeylSequencer
	/// Class for generating multivariate Weyl quasi-random sequences in arbitrary dimensions.
	struct WeylSequencer {

		// Base of the sequence
		vec<real> alphas;

		// Current state of the sequence
		vec<real> curr;


		/// Construct a multivariate Weyl sequence in the given dimensions.
		/// If no base is provided, the optimal choice is used.
		///
		/// @param dim Number of dimensions of the sequence
		/// @param alpha Base of the sequence
		WeylSequencer(size_t dim, real alpha = 0) {

			real base;
			alphas.resize(dim);

			if (alpha == 0) {

				// Find the only positive root of the polynomial: x^s+1 - x - 1 = 0
				base = 1.0 / root_bisect([=](real x) {
					return pow(x, dim + 1) - x - 1;
				}, 0.0, 2.0);

				// Optimal parameters are powers of the root inverse
				alpha = base;

			} else {
				base = alpha;
			}

			for (unsigned int i = 0; i < alphas.size(); ++i) {
				alphas[i] = alpha;
				alpha *= base;
			}

			curr = alphas;
		}


		/// Shift the Weyl sequence by the given amounts.
		template<typename Vector, enable_vector<Vector> = true>
		WeylSequencer& shift(const Vector& shifts) {

			if (alphas.size() != shifts.size()) {
				TH_MATH_ERROR("WeylSequencer::shift()", shifts.size(), MathError::InvalidArgument);
				return *this;
			}

			for (size_t i = 0; i < alphas.size(); i++)
				alphas[i] += shifts[i];

			return *this;
		}


		/// Overwrite a vector with the next point in the Weyl sequence.
		///
		/// @param v The vector to overwrite
		/// @return A reference to the modified vector
		template <typename Vector>
		Vector& operator()(Vector& v) {

			if (v.size() != alphas.size()) {
				TH_MATH_ERROR("WeylSequencer::operator()", v.size(), MathError::InvalidArgument);
				return algebra::vec_error(v);
			}

			for (size_t i = 0; i < v.size(); i++) {
				curr[i] = fract(curr[i] + alphas[i]);
				v[i] = curr[i];
			}
			
			return v;
		}


		/// Get the next point in the Weyl sequence.
		///
		/// @tparam The vector type to return
		/// @return A vector containing the next point
		template <typename Vector = vec<real>>
		Vector operator()() {

			Vector v;
			v.resize(alphas.size());
			return operator()(v);
		}
	};

}}

#endif
