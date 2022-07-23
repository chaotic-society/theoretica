
///
/// @file pseudorandom.h Pseudorandom number generation algorithms
///

#ifndef THEORETICA_PSEUDORANDOM_H
#define THEORETICA_PSEUDORANDOM_H

#include <vector>
#include <cstdint>
#include "../core/bit_op.h"

namespace theoretica {


    /// A pseudorandom function pointer
    using pseudorandom_function = uint64_t (*)(uint64_t, std::vector<uint64_t>&);


	/// Congruential pseudorandom number generation algorithm
	/// @param x The current recurrence value of the algorithm (x_n)
	/// @param a The multiplier term
	/// @param c The increment term
	/// @param m The modulus term
    /// @return The next generated pseudorandom number
	/// 
	/// The congruential generator is defined by the recurrence formula
	/// \f$x_{n+1} = (a x_n + c) mod m\f$ \n
    /// If no parameters are passed, the generator defaults to a = 48271,
    /// c = 0, m = (1 << 31) - 1.
	inline uint64_t rand_congruential(
        uint64_t x,
        uint64_t a = 48271,
        uint64_t c = 0,
        uint64_t m = ((uint64_t) 1 << 31) - 1) {

		return (a * x + c) % m;
	}


    /// Congruential pseudorandom number generation algorithm (wrapper)
    /// @param x The current recurrence value of the algorithm (x_n)
    /// @param state A vector containing the state of the algorithm (a, c, m in this order)
    /// @return The next generated pseudorandom number
    /// 
    /// \see rand_congruential
	inline uint64_t rand_congruential(uint64_t x, std::vector<uint64_t>& state) {

		if(state.size() != 3) {
			TH_MATH_ERROR("rand_congruential", state.size(), INVALID_ARGUMENT);
			return 0;
		}

		const uint64_t a = state[0];
		const uint64_t c = state[1];
		const uint64_t m = state[2];

        if(a > m || c > m) {
            TH_MATH_ERROR("rand_congruential", max(a, c), INVALID_ARGUMENT);
            return 0;
        }

		return rand_congruential(x, a, c, m);
	}


	/// Xoshiro256++ pseudorandom number generation algorithm
	/// @param x Dummy parameter (needed for function signature)
	/// @param state The four 64-bit integer state of the algorithm
	///
	/// Adapted from the reference implementation by Sebastiano Vigna
	inline uint64_t rand_xoshiro(uint64_t& a, uint64_t& b, uint64_t& c, uint64_t& d) {

		// Add and rotate
		const uint64_t result = bit_rotate(a + d, 23) + a;
		const uint64_t temp = b << 17;

		// Shift operations
		c ^= a;
		d ^= b;
		b ^= c;
		a ^= d;

		c ^= temp;
		d = bit_rotate(d, 45);

		return result;
	}


	/// Xoshiro256++ pseudorandom number generation algorithm
	/// @param x Dummy parameter (needed for function signature)
	/// @param state The four 64-bit integer state of the algorithm
	inline uint64_t rand_xoshiro(uint64_t x, std::vector<uint64_t>& state) {

		if(state.size() != 4) {
			TH_MATH_ERROR("rand_xoshiro", state.size(), INVALID_ARGUMENT);
			return 0;
		}

		return rand_xoshiro(state[0], state[1], state[2], state[3]);
	}


	/// SplitMix64 pseudorandom number generation
	/// @param x The 64-bit state of the algorithm
	///
	/// Adapted from the reference implementation by Sebastiano Vigna
	inline uint64_t rand_splitmix64(uint64_t x) {

		x += 0x9e3779b97f4a7c15;

		uint64_t res = x;
		res = (res ^ (res >> 30)) * 0xbf58476d1ce4e5b9;
		res = (res ^ (res >> 27)) * 0x94d049bb133111eb;

		return res ^ (res >> 31);
	}


	/// SplitMix64 pseudorandom number generation
	/// @param x The 64-bit state of the algorithm
	/// @param p Dummy variable (needed for function signature)
	inline uint64_t rand_splitmix64(uint64_t x, std::vector<uint64_t>& p) {
		return rand_splitmix64(x);
	}


	/// Wyrand pseudorandom number generation
	/// @param seed The (changing) seed of the algorithm
	/// @param p0 Additive constant (ideally a large prime number)
	/// @param p1 Mask for the algorithm
	///
	/// Algorithm by Yi Wang
	inline uint64_t rand_wyrand(uint64_t& seed, uint64_t p1, uint64_t p2) {
		seed += p1;
		return mix_mum(seed ^ p2, seed);
	}


	/// Wyrand pseudorandom number generation
	/// @param x Dummy variable
	/// @param p Algorithm parameters
	///
	/// p[0] is the initial seed, p[1] a large prime number and
	/// p[2] is the bit mask.
	inline uint64_t rand_wyrand(uint64_t x, std::vector<uint64_t>& p) {
		return rand_wyrand(p[0], p[1], p[2]);
	}

}


#endif
