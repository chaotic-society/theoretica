
///
/// @file pseudorandom_algo.h Pseudorandom number generation algorithms
///

#ifndef UROBORO_PSEUDO_RANDOM_H
#define UROBORO_PSEUDO_RANDOM_H

#include <vector>


namespace uroboro {


    /// A pseudorandom function pointer
    using pseudorandom_function = unsigned int (*)(unsigned int, const std::vector<unsigned int>&);


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
	unsigned int rand_congruential(
        unsigned int x,
        unsigned int a = 48271,
        unsigned int c = 0,
        unsigned int m = ((unsigned int) 1 << 31) - 1) {

		return (a * x + c) % m;
	}


    /// Congruential pseudorandom number generation algorithm (wrapper)
    /// @param x The current recurrence value of the algorithm (x_n)
    /// @param state A vector containing the state of the algorithm (a, c, m in this order)
    /// @return The next generated pseudorandom number
    /// 
    /// \see rand_congruential
	unsigned int rand_congruential(unsigned int x, const std::vector<unsigned int>& state) {

		if(state.size() != 3) {
			UMATH_ERROR("rand_congruential", state.size(), INVALID_ARGUMENT);
			return 0;
		}

		const unsigned int a = state[0];
		const unsigned int c = state[1];
		const unsigned int m = state[2];

        if(a > m || c > m) {
            UMATH_ERROR("rand_congruential", max(a, c), INVALID_ARGUMENT);
            return 0;
        }

		return rand_congruential(x, a, c, m);
	}

}


#endif
