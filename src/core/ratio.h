#ifndef UROBORO_RATIO_H
#define UROBORO_RATIO_H

namespace uroboro {

	template <typename T1, typename T2>
	class ratio {
		public:
			T1 num;
			T2 den;

			ratio() {
				num = 0;
				den = 0;
			}

			ratio(T1 n, T2 d){
				num = n;
				den = d;
			}

			~ratio() = default;

			ratio& operator*(ratio r) {
				return ratio(num * r.num, den * r.den);
			}

			ratio& operator/(ratio r) {
				return ratio(num * r.den, den * r.num);
			}

			ratio& operator+(ratio r) {
				return ratio(num * r.den + r.num * den, den * r.den);
			}

			ratio& operator-(ratio r) {
				return ratio(num * r.den - r.num * den, den * r.den);
			}

			ratio& operator*(T1 alpha) {
				return ratio(num * alpha, den);
			}

			ratio& operator/(T2 beta) {
				return ratio(num, den * beta);
			}


#ifndef UROBORO_NO_PRINT

			/// Convert the ratio to string representation
			inline std::string to_string() const {

				std::stringstream res;
				res << num << " / " << den;
				return res.str();
			}


			/// Stream the ratio in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const ratio& obj) {
				return out << obj.to_string();
			}

#endif


	};

}


#endif
