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
			std::vector<unsigned int> state;

		public:
			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and state s
			PRNG(pseudorandom_function p,
				unsigned int seed,
				const std::vector<unsigned int>& s) {

				f = p;
				x = seed;
				state = s;
			}

			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and state s
			PRNG(pseudorandom_function p,
				const std::vector<unsigned int>& s) {

				f = p;
				x = 1;
				state = s;
			}

			/// Seed the PRNG
			inline void seed(unsigned int seed) {
				x = seed;
			}


			/// Generate a pseudorandom number
			inline unsigned int next() {
				return x = f(x, state);
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


			/// Returns a standard linear congruential generator
			/// @param seed The seed to use for the generator (defaults to 1)
			/// @return A standard linear congruential PRNG object
			inline static PRNG linear_congruential(unsigned int seed = 1) {
				return PRNG(rand_congruential, seed, {48271, 0, ((unsigned int) 1 << 31) - 1});
			}
		
	};

}

#endif
