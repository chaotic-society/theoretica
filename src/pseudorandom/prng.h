
///
/// @file prng.h Pseudorandom number generator class
///

#ifndef UROBORO_PRNG_H
#define UROBORO_PRNG_H

#include <vector>


namespace uroboro {

	/// @class PRNG
	/// A pseudorandom number generator
	class PRNG {
		private:
			pseudorandom_function f;
			unsigned int x;
			std::vector<unsigned int> param;

		public:
			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and parameters s
			PRNG(pseudorandom_function p,
				unsigned int seed,
				const std::vector<unsigned int>& s) {

				f = p;
				x = seed;
				param = s;
			}

			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and parameters s
			PRNG(pseudorandom_function p,
				const std::vector<unsigned int>& s) {

				f = p;
				x = 1;
				param = s;
			}

			/// Seed the PRNG
			inline void seed(unsigned int seed) {
				x = seed;
			}


			/// Generate a pseudorandom number
			inline unsigned int next() {
				return x = f(x, param);
			}


			/// Generate a pseudorandom number
			/// \see next()
			inline unsigned int operator()() {
				return next();
			}


			/// Discard n numbers from the generator.
			/// Equivalent to calling next() n times.
			inline void discard(unsigned int n) {
				for (int i = 0; i < n; ++i)
					next();
			}


			/// Set the generating function
			inline void set_function(pseudorandom_function p) {
				f = p;
			}

			/// Get the generating function
			inline pseudorandom_function get_function() {
				return f;
			}


			/// Set the generator's parameters
			inline void set_param(std::vector<unsigned int> v) {
				param = v;
			}

			/// Get the generator's parameters
			inline std::vector<unsigned int> get_param() {
				return param;
			}


			/// Returns a standard linear congruential generator
			/// @param seed The seed to use for the generator (defaults to 1)
			/// @return A standard linear congruential PRNG object
			inline static PRNG linear_congruential(unsigned int seed = 1) {
				return PRNG(rand_congruential, seed, {48271, 0, ((unsigned int) 1 << 31) - 1});
			}
		
	};


	/// Generate a pseudorandom real number in [a, b] using a
	/// preexisting generator.
	///
	/// The algorithm generates natural numbers between 0 and prec,
	/// scales the result to [0, 1] and then transforms it to [a, b]
	inline real rand_real(real a, real b, PRNG& g, unsigned int prec = (1 << 30)) {
		return a + (b - a) * (g.next() % prec) / static_cast<real>(prec);
	}

}

#endif
