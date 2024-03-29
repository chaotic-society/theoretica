
///
/// @file polynomial.h Polynomial storage and manipulation
///

#ifndef THEORETICA_POLYNOMIAL_H
#define THEORETICA_POLYNOMIAL_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/vec_buff.h"
#include "../core/real_analysis.h"
#include "../algebra/vec.h"
#include "../complex/complex.h"
#include "../complex/complex_analysis.h"


namespace theoretica {

	/// @class polynomial
	/// A polynomial of arbitrary order
	template<typename T = real>
	class polynomial {
		public:
			std::vector<T> coeff;

			/// Initialize as an empty polynomial
			polynomial() : coeff() {}

			/// Initialize as a constant
			polynomial(T a) : coeff({a}) {}

			/// Initialize from an std::vector
			polynomial(const std::vector<T>& c) : coeff(c) {}
					
			/// Initialize from an std::initializer_list
			polynomial(std::initializer_list<T> l) : coeff(l) {}

			/// Default destructor
			~polynomial() {}


			/// Access i-th coefficient
			inline T& at(int i) {
				return coeff[i];
			}


			/// Get the i-th by value
			inline T get(int i) const {
				return coeff[i];
			}


			/// Return the nth order coefficient
			inline T& operator[](unsigned int i) {
				return coeff[i];
			}


			/// Evaluate the polynomial using x as variable
			template<typename EvalType = T>
			inline EvalType eval(EvalType x) const {

				if(!coeff.size())
					return 0;

				EvalType sum = 0;

				// Evaluate using Horner's method
				for (unsigned int i = 0; i < coeff.size(); ++i)
					sum = coeff[coeff.size() - i - 1] + x * sum;

				// TO-DO Compare performance using fma x86 Assembly instruction
				// sum = fma(x, sum, coeff[coeff.size() - i - 1]);

				return sum;
			}


			/// Evaluate the polynomial using x as variable
			template<typename EvalType = T>
			inline EvalType operator()(EvalType x) const {
				return eval(x);
			}


			/// Find the true order of the polynomial (ignoring null coefficients)
			inline unsigned int find_order() const {

				for (int i = coeff.size() - 1; i >= 0; --i) {
					if(coeff[i] != 0)
						return i;
				}

				return 0;
			}


			/// Remove higher order null coefficients
			inline void trim() {

				for (unsigned int i = coeff.size() - 1; i >= 0; --i) {
					if(abs(coeff[i]) > MACH_EPSILON)
						break;

					coeff.pop_back();
				}
			}


			/// Get the number of coefficients
			inline size_t size() const {
				return coeff.size();
			}


			/// Sum two polynomials
			inline polynomial operator+(const polynomial& p) const {

				polynomial r(coeff);

				for (unsigned int i = 0; i < min(r.size(), p.size()); ++i) {
					r[i] += p.get(i);
				}

				return r;
			}


			/// Subtract a polynomial from another
			inline polynomial operator-(const polynomial& p) const {
				
				polynomial r(coeff);

				for (unsigned int i = 0; i < min(r.size(), p.size()); ++i) {
					r[i] -= p.get(i);
				}

				return r;
			}


			/// Multiply two polynomials
			inline polynomial operator*(const polynomial& p) const {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1);

				for (unsigned int i = 0; i < size(); ++i) {
					for (unsigned int j = 0; j < p.size(); ++j) {
						r[i + j] += coeff[i] * p.get(j);
					}
				}

				return r;
			}


			/// Polynomial division
			inline polynomial operator/(const polynomial& d) const {

				const unsigned int d_order = d.find_order();
				const unsigned int this_order = find_order();

				if(d_order == 0 && d.get(0) == 0) {
					TH_MATH_ERROR("polynomial::operator/", d.get(0), DIV_BY_ZERO);
					return polynomial(nan());
				}

				// Remainder
				polynomial r = *this;

				// Quotient
				polynomial q = polynomial();
				unsigned int i = 0;

				while(i < this_order) {

					// Compute only once the degree of the polynomial
					const unsigned int r_order = r.find_order();

					// Stop execution if the division is complete
					// (when the remainder is 0 or has lower degree)
					if((r_order == 0 && (abs(r.get(0)) < MACH_EPSILON)) || r_order < d_order)
						break;

					// Simple division between highest degree terms
					const polynomial t = polynomial<T>::monomial(
						r.get(r_order) / d.get(d_order),
						r_order - d_order);

					// Add monomial to quotient and subtract the
					// monomial times the dividend from the remainder
					q += t;
					r -= t * d;

					i++;
				}

				// The algorithm has stopped iterating after a number
				// of the dividend's degree counts
				if(i == this_order) {
					TH_MATH_ERROR("polynomial::operator/", i, NO_ALGO_CONVERGENCE);
					return polynomial(nan());
				}

				return q;
			}


			/// Multiply a polynomial by a scalar
			inline polynomial operator*(T a) const {

				polynomial r = polynomial(*this);

				for (unsigned int i = 0; i < size(); ++i)
					r.coeff[i] *= a;

				return r;
			}


			/// Divide a polynomial by a scalar
			inline polynomial operator/(T a) const {

				if(a == 0) {
					TH_MATH_ERROR("polynomial::operator/", a, DIV_BY_ZERO);
				}

				polynomial r = polynomial(*this);

				for (unsigned int i = 0; i < size(); ++i)
					r.coeff[i] /= a;

				return r;
			}


			/// Sum a polynomial to this one
			inline polynomial& operator+=(const polynomial& p) {

				// Make room for the new coefficients
				if(coeff.size() < p.size())
					coeff.resize(p.size(), T(0));

				for (unsigned int i = 0; i < min(size(), p.size()); ++i)
					coeff[i] += p.get(i);

				return *this;
			}


			/// Subtract a polynomial from this one
			inline polynomial& operator-=(const polynomial& p) {

				// Make room for the new coefficients
				if(coeff.size() < p.size())
					coeff.resize(p.size(), T(0));

				for (unsigned int i = 0; i < min(size(), p.size()); ++i)
					coeff[i] -= p.get(i);

				return *this;
			}


			/// Multiply two polynomials
			inline polynomial& operator*=(const polynomial& p) {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1, T(0));

				for (unsigned int i = 0; i < size(); ++i)
					for (unsigned int j = 0; j < p.size(); ++j)
						r[i + j] += coeff[i] * p.get(j);

				*this = r;
				return *this;
			}


			/// Multiply a polynomial by a scalar value
			inline polynomial& operator*=(T a) {

				for (unsigned int i = 0; i < coeff.size(); ++i)
					coeff[i] *= a;
				
				return *this;
			}


			/// Multiply a polynomial by a scalar value
			inline polynomial& operator/=(const polynomial& a) {
				return (*this = (*this / a));
			}


			/// Divide a polynomial by a scalar value
			inline polynomial& operator/=(T a) {

				if(a == 0) {
					TH_MATH_ERROR("polynomial::operator/=", a, DIV_BY_ZERO);
					return *this;
				}

				for (unsigned int i = 0; i < coeff.size(); ++i)
					coeff[i] /= a;
				
				return *this;
			}


			/// Check whether two polynomials are equal
			inline bool operator==(const polynomial& other) const {

				for (unsigned int i = 0; i < min(other.size(), this->size()); ++i) {
					if(other.coeff[i] != coeff[i])
						return false;
				}

				if(size() > other.size()) {
					for (unsigned int i = 0; i < size(); ++i) {
						if(coeff[i] != 0)
							return false;
					}
				} else {
					for (unsigned int i = 0; i < other.size(); ++i) {
						if(other.coeff[i] != 0)
							return false;
					}
				}

				return true;
			}


			/// Compute the roots of a quadratic polynomial
			inline vec<complex<>, 2> quadratic_roots() const {

				const int order = find_order();

				// Check that the polynomial is quadratic
				if(order != 2) {
					TH_MATH_ERROR("quadratic_roots", order, IMPOSSIBLE_OPERATION);
					return vec<complex<>, 2>({nan(), nan()});
				}

				const T p = coeff[1] / coeff[2];
				const T q = coeff[0] / coeff[2];

				// Case when 0 is a root
				if(abs(q) < MACH_EPSILON)
					return {complex<>(-p), complex<>(0)};

				complex<> z1, z2;

				// Use Vieta's theorem to avoid catastrophic cancellation
				if(abs(p) > 1) {

					z1 = -sgn(p) * (abs(p) / 2.0
						+ abs(p) * sqrt(complex<>(0.25 - (q / p) / p)));
					z2 = q / z1;

				} else {

					const complex<> s = sqrt(complex<>(0.25 * square(p) - q));
					z1 = -p / 2.0 + s;
					z2 = -p / 2.0 - s;
				}

				return {z1, z2};
			}


			/// Construct a polynomial from its roots
			inline static polynomial<T> from_roots(const std::vector<T>& roots) {

				polynomial<T> P = {1};

				// P = product((x - r_i))
				for (unsigned int i = 0; i < roots.size(); ++i)
					P *= polynomial<T>({roots[i] * -1, 1});

				return P;
			}


			/// Returns a monomial of the given degree and coefficient
			inline static polynomial<T> monomial(T c, unsigned int order) {
				
				polynomial m;
				m.coeff = std::vector<T>(order + 1, T(0));
				m.coeff[order] = c;
				return m;
			}


			// Friend operators to enable equations of the form
			// (T) op. (polynomial<T>)

			inline friend polynomial<T> operator+(T r, const polynomial<T>& z) {
				return z + polynomial(r);
			}

			inline friend polynomial<T> operator-(T r, const polynomial<T>& z) {
				return (z * -1) + polynomial(r);
			}

			inline friend polynomial<T> operator*(T r, const polynomial<T>& z) {
				return z * r;
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the polynomial to string representation
			inline std::string to_string(
				const std::string& unknown = "x",
				const std::string& exponentiation = "^") const {

				std::stringstream res;
				bool flag = false;
				const int sz = coeff.size();

				for (int i = sz - 1; i >= 0; --i) {

					if(abs(coeff[i]) < MACH_EPSILON)
						continue;

					res << (coeff[i] >= 0 ? "+ " : "- ");
					res << abs(coeff[i]);

					if(i) {
						res << "*" << unknown << exponentiation << i;
						res << " ";
					}

					flag = true;
				}

				if(!flag)
					res << "0";

				return res.str();
			}


			/// Stream the polynomial in string representation
			/// to an output stream (std::ostream)
			friend std::ostream& operator<<(std::ostream& out, const polynomial& obj) {
				return out << obj.to_string();
			}

#endif

		
	};
	
}

#endif
