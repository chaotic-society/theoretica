
///
/// @file prng.h Pseudorandom number generation.
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
			
			/// A function which takes the state of the generator
			/// and returns the next generated pseudorandom number.
			pseudorandom_function f;

			/// The last generated pseudorandom number.
			uint64_t x;

			/// The state of the pseudorandom generator.
			std::vector<uint64_t> param;

		public:

			/// Construct a PRNG with the given
			/// generating algorithm p, seed x and parameters s
			PRNG(pseudorandom_function p,
				uint64_t seed,
				const std::vector<uint64_t>& s) : f(p), x(seed), param(s) {}


			/// Construct a PRNG with the given
			/// generating algorithm and seed
			PRNG(pseudorandom_function p, uint64_t seed) : f(p), x(seed) {}


			/// Construct a PRNG with the default
			/// generator and the given seed.
			PRNG(uint64_t seed) {
				*this = PRNG::xoshiro(seed);
			}


			/// Seed the PRNG.
			inline void seed(uint64_t seed) {
				x = seed;
			}


			/// Generate a pseudorandom number.
			inline uint64_t next() {
				return x = f(x, param);
			}


			/// Generate a pseudorandom number.
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


			/// Return the last generated number.
			inline uint64_t last() const {
				return x;
			}


			/// Set the generating function.
			inline void set_function(pseudorandom_function p) {
				f = p;
			}

			/// Get the generating function.
			inline pseudorandom_function get_function() const {
				return f;
			}


			/// Set the generator's parameters.
			inline void set_param(const std::vector<uint64_t>& v) {
				param = v;
			}

			/// Set a specific parameter by index.
			inline void set_param(unsigned int i, uint64_t value) {
				param[i] = value;
			}

			/// Get the generator's parameters.
			inline std::vector<uint64_t> get_param() const {
				return param;
			}


			/// Stream the next generated number
			inline PRNG& operator>>(uint64_t& n) {
				n = next();
				return *this;
			}


			/// Returns a standard linear congruential generator
			/// @param seed The seed to use for the generator (defaults to 1)
			/// @return A standard linear congruential PRNG object
			inline static PRNG linear_congruential(uint64_t seed = 1) {

				if(seed == 0)
					seed = 1;

				return PRNG(
					randgen_congruential, seed,
					{48271, 0, ((uint64_t) 1 << 31) - 1}
				);
			}


			/// Returns a Xoshiro256++ generator
			/// @param p The four state parameters
			/// @return A Xoshiro256++ pseudorandom number generator
			inline static PRNG xoshiro(const std::vector<uint64_t>& p) {
				return PRNG(randgen_xoshiro, 0, p);
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
				
				const uint64_t n1 = randgen_splitmix64(seed);
				const uint64_t n2 = randgen_splitmix64(n1);
				const uint64_t n3 = randgen_splitmix64(n2);
				const uint64_t n4 = randgen_splitmix64(n3);

				return PRNG(randgen_xoshiro, 0, {n1, n2, n3, n4});
			}


			/// Returns a Splitmix64 generator
			/// @param seed The seed to use for the generator
			/// @return A splitmix64 pseudorandom number generator
			inline static PRNG splitmix64(uint64_t seed = 1) {

				if(seed == 0)
					seed = 1;

				return PRNG(randgen_splitmix64, seed);
			}


			/// Returns a Wyrand generator
			inline static PRNG wyrand(uint64_t seed = 1,
				uint64_t p1 = 2549536629329, uint64_t p2 = 136137137) {

				if(seed == 0)
					seed = 1;

				if(p2 == 0)
					p2 = randgen_splitmix64(seed);

				return PRNG(randgen_wyrand, 0, {seed, p1, p2});
			}


			/// Returns a Middle-square generator
			inline static PRNG middlesquare(
				uint64_t seed, uint64_t offset = 765872292751861) {

				if(seed == 0)
					seed = randgen_splitmix64(765872292751861);

				if(offset == 0)
					offset = 765872292751861;

				return PRNG(randgen_middlesquare, seed, {offset});
			}
		
	};


	/// Shuffle a set by exchanging random couples of elements.
	/// 
	/// @param v The set to shuffle (as a Vector with [] and .size() methods)
	/// @param g An already initialized PRNG
	/// @param rounds The number of pairs to exchange.
	/// @return A reference to the shuffled vector
	template<typename Vector>
	inline Vector& shuffle(Vector& v, PRNG& g, unsigned int rounds) {

		if(v.size() == 0)
			TH_MATH_ERROR("shuffle", v.size(), INVALID_ARGUMENT);

		for (unsigned int i = 0; i < rounds; ++i) {
			
			// Generate two random indices
			const size_t index1 = g() % v.size();
			const size_t index2 = g() % v.size();

			// Exchange the two values at the random indices
			const auto x1 = v[index1];
			v[index1] = v[index2];
			v[index2] = x1;
		}
	}


	/// Shuffle a set by exchanging random couples of elements.
	/// The number of rounds used is \f$(N - 1)^2\f$.
	/// 
	/// @param v The set to shuffle (as a Vector with [] and .size() methods)
	/// @param g An already initialized PRNG
	/// @return A reference to the shuffled vector
	template<typename Vector>
	inline Vector& shuffle(Vector& v, PRNG& g) {
		return shuffle(v, g, square(v.size() - 1));
	}

}

#endif
