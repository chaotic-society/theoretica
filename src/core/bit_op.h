
///
/// @file bit_op.h Operations on bits
///

#ifndef THEORETICA_BIT_OP_H
#define THEORETICA_BIT_OP_H

#include <cstdint>


namespace theoretica {


	/// Multiply two 64-bit integers and store the result
	/// in two 64-bit variables, keeping 128 bits of the result.
	/// @param a The first number to multiply
	/// @param b The second number to multiply
	/// @param c_low The variable where to store the lowest 64 bits
	/// of the result
	/// @param c_high The variable where to store the highest 64 bits
	/// of the result
	inline void mul_uint128(
		uint64_t a, uint64_t b,
		uint64_t& c_low, uint64_t& c_high) {

		// Lowest and highest 32 bits of a and b
		uint64_t a_low = a & 0xffffffff;
		uint64_t a_high = a >> 32;

		uint64_t b_low = b & 0xffffffff;
		uint64_t b_high = b >> 32;

		uint64_t m[4];

		// Multiplication terms for (a_l + a_h) * (b_l + b_h)
		m[0] = a_low * b_low;
		m[1] = a_low * b_high;
		m[2] = a_high * b_low;
		m[3] = a_high * b_high;

		// Multiplication carry for c_high
		uint64_t carry = ((m[0] >> 32) + (m[1] & 0xffffffff) + (m[2] & 0xffffffff)) >> 32;

		c_low = m[0] + (m[1] << 32) + (m[2] << 32);
		c_high = m[3] + (m[1] >> 32) + (m[2] >> 32) + carry;
	}


	/// MUM bit mixing function
	/// @param a The first operand
	/// @param b The second operand
	inline uint64_t mix_mum(uint64_t a, uint64_t b) {

		uint64_t c_low, c_high;
		mul_uint128(a, b, c_low, c_high);

		return c_high ^ c_low;
	}


	/// Bit rotation using shifting
	/// @param x The variable to rotate the bits of
	/// @param i The index of the rotated bits
	template<typename ShiftableType>
	inline ShiftableType bit_rotate(ShiftableType x, unsigned int i) {
		return (x << i) | (x >> ((sizeof(ShiftableType) * 8) - i));
	}

}


#endif
