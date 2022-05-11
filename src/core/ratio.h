
///
/// @file ratio.h A ratio between objects
///

#ifndef UROBORO_RATIO_H
#define UROBORO_RATIO_H

namespace uroboro {

	/// @class ratio A class representing a ratio between two objects,
	/// like a fraction or a rational polynomial.
	///
	/// @note For the class to properly work, the template types
	/// need to have operator*, operator+ and operator- defined.
	template <typename T1, typename T2>
	class ratio {
		public:

			/// The numerator
			T1 num;

			/// The denominator
			T2 den;

			/// Construct the object from the default contructors
			/// of the two types
			ratio() {
				num = T1();
				den = T2();
			}

			/// Construct the object from two objects
			///
			/// @param n The numerator
			/// @param d The denominator
			ratio(T1 n, T2 d){
				num = n;
				den = d;
			}

			~ratio() = default;

			/// Multiply two ratios
			inline ratio operator*(ratio r) const {
				return ratio(num * r.num, den * r.den);
			}

			/// Divide two ratios (without explicitly using division)
			inline ratio operator/(ratio r) const {
				return ratio(num * r.den, den * r.num);
			}

			/// Add two ratios
			inline ratio operator+(ratio r) const {
				return ratio(num * r.den + r.num * den, den * r.den);
			}

			/// Subtract two ratios
			inline ratio operator-(ratio r) const {
				return ratio(num * r.den - r.num * den, den * r.den);
			}

			/// Multiply the ratio by a factor
			inline ratio operator*(T1 a) const {
				return ratio(num * a, den);
			}

			/// Divide the ratio by a factor
			inline ratio operator/(T2 b) const {
				return ratio(num, den * b);
			}

			/// Evaluate the ratio as the division
			/// between numerator and denominator converted
			/// to the specified type.
			///
			/// This function is well-defined only if division
			/// between objects of type T is defined.
			/// A static cast is performed before division.
			template<typename T>
			inline T eval_as() {
				return static_cast<T>(num) / static_cast<T>(den);
			}


#ifndef UROBORO_NO_PRINT

			/// Convert the ratio to string representation
			inline std::string to_string() const {

				std::stringstream res;
				res << num << "/" << den;
				return res.str();
			}


			/// Stream the ratio in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const ratio<T1, T2>& obj) {
				return out << obj.to_string();
			}

#endif

	};

}


#endif
