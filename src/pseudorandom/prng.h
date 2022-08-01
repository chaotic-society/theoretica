
///
/// @file prng.h Pseudorandom number generator class
///

#ifndef THEORETICA_PRNG_H
#define THEORETICA_PRNG_H

#include <vector>
#include "./pseudorandom.h"
#include "../core/error.h"


namespace theoretica {

	/// @class PRNG
	/// A pseudorandom number generator
	class PRNG {
		private:
			pseudorandom_function f;
			uint64_t x;
			std::vector<uint64_t> param;

		public:
			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and parameters s
			PRNG(pseudorandom_function p,
				uint64_t seed,
				const std::vector<uint64_t>& s) : f(p), x(seed), param(s) {}

			/// Construct a PRNG with the given
			/// generating algorithm p and parameters s
			///
			/// The seed will be set to 1.
			PRNG(pseudorandom_function p,
				const std::vector<uint64_t>& s) : f(p), x(1), param(s) {}

			/// Construct a PRNG with the given
			/// generating algorithm and seed
			PRNG(pseudorandom_function p, uint64_t seed) : f(p), x(seed) {}

			/// Seed the PRNG
			inline void seed(uint64_t seed) {
				x = seed;
			}


			/// Generate a pseudorandom number
			inline uint64_t next() {
				return x = f(x, param);
			}


			/// Generate a pseudorandom number
			/// \see next()
			inline uint64_t operator()() {
				return next();
			}


			/// Discard n numbers from the generator.
			/// Equivalent to calling next() n times.
			inline void discard(uint64_t n) {
				for (uint64_t i = 0; i < n; ++i)
					next();
			}


			/// Return the last generated number
			inline uint64_t last() const {
				return x;
			}


			/// Set the generating function
			inline void set_function(pseudorandom_function p) {
				f = p;
			}

			/// Get the generating function
			inline pseudorandom_function get_function() const {
				return f;
			}


			/// Set the generator's parameters
			inline void set_param(const std::vector<uint64_t>& v) {
				param = v;
			}

			/// Set a specific parameter
			inline void set_param(unsigned int i, uint64_t value) {
				param[i] = value;
			}

			/// Get the generator's parameters
			inline std::vector<uint64_t> get_param() const {
				return param;
			}


			/// Returns a standard linear congruential generator
			/// @param seed The seed to use for the generator (defaults to 1)
			/// @return A standard linear congruential PRNG object
			inline static PRNG linear_congruential(uint64_t seed = 1) {

				if(seed == 0)
					seed = 1;

				return PRNG(rand_congruential, seed, {48271, 0, ((uint64_t) 1 << 31) - 1});
			}


			/// Returns a Xoshiro256++ generator
			/// @param p The four state parameters
			/// @return A Xoshiro256++ pseudorandom number generator
			inline static PRNG xoshiro(const std::vector<uint64_t>& p) {
				return PRNG(rand_xoshiro, 0, p);
			}


			/// Returns a Xoshiro256++ generator
			/// @param seed The seed to use for the generator
			/// @return A Xoshiro256++ pseudorandom number generator
			///
			/// The four parameters for the Xoshiro256++ algorithm
			/// are generated using the Splitmix64 algorithm.
			inline static PRNG xoshiro(uint64_t seed = 1) {

				if(seed == 0)
					seed = 1;
				
				const uint64_t n1 = rand_splitmix64(seed);
				const uint64_t n2 = rand_splitmix64(n1);
				const uint64_t n3 = rand_splitmix64(n2);
				const uint64_t n4 = rand_splitmix64(n3);

				return PRNG(rand_xoshiro, 0, {n1, n2, n3, n4});
			}


			/// Returns a Splitmix64 generator
			/// @param seed The seed to use for the generator
			/// @return A splitmix64 pseudorandom number generator
			inline static PRNG splitmix64(uint64_t seed = 1) {

				if(seed == 0)
					seed = 1;

				return PRNG(rand_splitmix64, seed);
			}


			/// Returns a Wyrand generator
			inline static PRNG wyrand(uint64_t seed = 1,
				uint64_t p1 = 2549536629329, uint64_t p2 = 136137137) {

				if(seed == 0)
					seed = 1;

				if(p2 == 0)
					p2 = rand_splitmix64(seed);

				return PRNG(rand_wyrand, 0, {seed, p1, p2});
			}
		
	};

}

#endif
