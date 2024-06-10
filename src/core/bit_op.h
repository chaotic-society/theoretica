
///
/// @file bit_op.h Operations on bits
///

#ifndef THEORETICA_BIT_OP_H
#define THEORETICA_BIT_OP_H

#include <cstdint>


namespace theoretica {


	/// Multiply two 64-bit unsigned integers and store the result
	/// in two 64-bit variables, keeping 128 bits of the result.
	///
	/// @param a The first number to multiply
	/// @param b The second number to multiply
	/// @param c_low The variable to store the lowest 64 bits
	/// of the result.
	/// @param c_high The variable to store the highest 64 bits
	/// of the result.
	inline void mul_uint128(
		uint64_t a, uint64_t b,
		uint64_t& c_low, uint64_t& c_high) {

		// Lowest and highest 32 bits of a and b
		const uint64_t a_low = a & 0xffffffff;
		const uint64_t a_high = a >> 32;

		const uint64_t b_low = b & 0xffffffff;
		const uint64_t b_high = b >> 32;

		uint64_t m[4];

		// Multiplication terms for (a_l + a_h) * (b_l + b_h)
		m[0] = a_low * b_low;
		m[1] = a_low * b_high;
		m[2] = a_high * b_low;
		m[3] = a_high * b_high;

		// Multiplication carry for c_high
		uint64_t carry = (
			(m[0] >> 32) +
			(m[1] & 0xffffffff) +
			(m[2] & 0xffffffff)) >> 32;

		c_low = m[0] + (m[1] << 32) + (m[2] << 32);
		c_high = m[3] + (m[1] >> 32) + (m[2] >> 32) + carry;
	}


	/// MUM bit mixing function, computes the 128-bit
	/// product of a and b and the XOR of their high
	/// and low 64-bit parts.
	///
	/// @param a The first operand
	/// @param b The second operand
	/// @return The XOR of the high and low bits
	/// of the 128-bit product of a and b.
	inline uint64_t mix_mum(uint64_t a, uint64_t b) {

		uint64_t c_low, c_high;
		mul_uint128(a, b, c_low, c_high);

		return c_high ^ c_low;
	}


	/// Bit rotation of unsigned integer types using shifts.
	///
	/// @param x The unsigned integer to rotate the bits of
	/// @param i The index of the rotated bits
	/// @return The unsigned integer with the given bits rotated
	template<typename UnsignedIntType>
	inline UnsignedIntType bit_rotate(UnsignedIntType x, unsigned int i) {

		return (x << i) | (x >> ((sizeof(UnsignedIntType) * 8) - i));
	}

}


#endif
