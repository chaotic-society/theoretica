
///
/// @file polynomial.h Polynomial storage and manipulation
///

#ifndef THEORETICA_POLYNOMIAL_H
#define THEORETICA_POLYNOMIAL_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/real_analysis.h"
#include "../algebra/vec.h"
#include "../complex/complex.h"
#include "../complex/complex_analysis.h"


namespace theoretica {

	/// @class polynomial
	/// A polynomial of arbitrary order
	template<typename Type = real>
	class polynomial {
		public:
			std::vector<Type> coeff;

			/// Initialize as an empty polynomial
			polynomial() : coeff() {}

			/// Initialize as a constant
			polynomial(Type a) : coeff({a}) {}

			/// Initialize from an std::vector
			polynomial(const std::vector<Type>& c) : coeff(c) {}
					
			/// Initialize from an std::initializer_list
			polynomial(std::initializer_list<Type> l) : coeff(l) {}

			/// Default destructor
			~polynomial() {}


			/// Get i-th coefficient by constant reference, with bound checking.
			inline const Type& at(unsigned int i) const {
				return coeff.at(i);
			}

			/// Access i-th coefficient by reference, with bound checking.
			inline Type& at(unsigned int i) {
				return coeff.at(i);
			}


			/// Get the n-th order coefficient by constant reference.
			inline const Type& operator[](unsigned int i) const {
				return coeff[i];
			}


			/// Get the n-th order coefficient by reference.
			inline Type& operator[](unsigned int i) {
				return coeff[i];
			}


			/// Evaluate the polynomial using x as variable
			template<typename EvalType = Type>
			inline EvalType eval(EvalType x) const {

				EvalType sum = EvalType(0.0);

				// Evaluate using Horner's method
				for (unsigned int i = 0; i < coeff.size(); ++i)
					sum = coeff[coeff.size() - i - 1] + x * sum;

				return sum;
			}


			/// Evaluate the polynomial using x as variable
			template<typename EvalType = Type>
			inline EvalType operator()(EvalType x) const {
				return eval(x);
			}


			/// Find the true order of the polynomial
			/// (ignoring trailing null coefficients)
			inline unsigned int find_order() const {

				for (int i = coeff.size() - 1; i >= 0; --i) {
					if(abs(coeff[i]) < MACH_EPSILON)
						return i;
				}

				return 0;
			}


			/// Remove trailing zero coefficients
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

				polynomial r;
				r.coeff.resize(max(coeff.size(), p.size()));

				unsigned int i = 0;

				for (; i < min(coeff.size(), p.size()); ++i)
					r[i] = coeff[i] + p[i];

				if (coeff.size() > p.size())
					for (; i < coeff.size(); ++i)
						r[i] = coeff[i];
				else
					for (; i < p.size(); ++i)
						r[i] = p[i];

				return r;
			}


			/// Subtract a polynomial from another
			inline polynomial operator-(const polynomial& p) const {
				
				polynomial r;
				r.coeff.resize(max(coeff.size(), p.size()));

				unsigned int i = 0;

				for (; i < min(r.size(), p.size()); ++i)
					r[i] = coeff[i] - p[i];

				if (coeff.size() > p.size())
					for (; i < coeff.size(); ++i)
						r[i] = coeff[i];
				else
					for (; i < p.size(); ++i)
						r[i] = -p[i];

				return r;
			}


			/// Multiply two polynomials
			inline polynomial operator*(const polynomial& p) const {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1);

				for (unsigned int i = 0; i < size(); ++i) {
					for (unsigned int j = 0; j < p.size(); ++j) {
						r[i + j] += coeff[i] * p[j];
					}
				}

				return r;
			}


			/// Polynomial division
			inline polynomial operator/(const polynomial& d) const {

				const unsigned int d_order = d.find_order();
				const unsigned int this_order = find_order();

				if(d_order == 0 && d[0] == 0) {
					TH_MATH_ERROR("polynomial::operator/", d[0], DIV_BY_ZERO);
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
					if((r_order == 0 && (abs(r[0]) < MACH_EPSILON)) || r_order < d_order)
						break;

					// Simple division between highest degree terms
					const polynomial t = polynomial<Type>::monomial(
						r[r_order] / d[d_order],
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
			inline polynomial operator*(Type a) const {

				polynomial r = polynomial(*this);

				for (unsigned int i = 0; i < size(); ++i)
					r.coeff[i] *= a;

				return r;
			}


			/// Divide a polynomial by a scalar
			inline polynomial operator/(Type a) const {

				if(abs(a) < MACH_EPSILON) {
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
					coeff.resize(p.size(), Type(0));

				for (unsigned int i = 0; i < min(size(), p.size()); ++i)
					coeff[i] += p[i];

				return *this;
			}


			/// Subtract a polynomial from this one
			inline polynomial& operator-=(const polynomial& p) {

				// Make room for the new coefficients
				if(coeff.size() < p.size())
					coeff.resize(p.size(), Type(0));

				for (unsigned int i = 0; i < min(coeff.size(), p.size()); ++i)
					coeff[i] -= p[i];

				return *this;
			}


			/// Multiply two polynomials
			inline polynomial& operator*=(const polynomial& p) {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1, Type(0));

				for (unsigned int i = 0; i < coeff.size(); ++i)
					for (unsigned int j = 0; j < p.size(); ++j)
						r[i + j] += coeff[i] * p[j];

				*this = r;
				return *this;
			}


			/// Multiply a polynomial by a scalar value
			inline polynomial& operator*=(Type a) {

				for (unsigned int i = 0; i < coeff.size(); ++i)
					coeff[i] *= a;
				
				return *this;
			}


			/// Multiply a polynomial by a scalar value
			inline polynomial& operator/=(const polynomial& a) {
				return (*this = (*this / a));
			}


			/// Divide a polynomial by a scalar value
			inline polynomial& operator/=(Type a) {

				if(abs(a) < MACH_EPSILON) {
					TH_MATH_ERROR("polynomial::operator/=", a, DIV_BY_ZERO);
					return *this;
				}

				for (unsigned int i = 0; i < coeff.size(); ++i)
					coeff[i] /= a;
				
				return *this;
			}


			/// Check whether two polynomials are equal
			inline bool operator==(const polynomial& other) const {

				const unsigned int n = min(other.size(), this->size());
				for (unsigned int i = 0; i < n; ++i) {
					if(abs(other.coeff[i] - coeff[i]) >= MACH_EPSILON)
						return false;
				}

				if(size() > other.size()) {
					for (unsigned int i = 0; i < size(); ++i) {
						if(abs(coeff[i]) >= MACH_EPSILON)
							return false;
					}
				} else {
					for (unsigned int i = 0; i < other.size(); ++i) {
						if(abs(other.coeff[i]) >= MACH_EPSILON)
							return false;
					}
				}

				return true;
			}


			/// Get an iterator for the first coefficient
			/// of the polynomial.
			inline auto begin() {
				return coeff.begin();
			}


			/// Get an iterator for the one plus last coefficient
			/// of the polynomial.
			inline auto end() {
				return coeff.end();
			}


			/// Get a const iterator for the first coefficient
			/// of the polynomial.
			inline auto cbegin() {
				return coeff.cbegin();
			}


			/// Get a const iterator for the one plus last coefficient
			/// of the polynomial.
			inline auto cend() {
				return coeff.cend();
			}


			/// Compute the roots of a quadratic polynomial
			inline vec<complex<>, 2> quadratic_roots() const {

				const int order = find_order();

				// Check that the polynomial is quadratic
				if(order != 2) {
					TH_MATH_ERROR("quadratic_roots", order, IMPOSSIBLE_OPERATION);
					return vec<complex<>, 2>({nan(), nan()});
				}

				const Type p = coeff[1] / coeff[2];
				const Type q = coeff[0] / coeff[2];

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
			inline static polynomial<Type> from_roots(
				const std::vector<Type>& roots) {

				polynomial<Type> P = {1};

				// P = product((x - r_i))
				for (unsigned int i = 0; i < roots.size(); ++i)
					P *= polynomial<Type>({roots[i] * -1, 1});

				return P;
			}


			/// Returns a monomial of the given degree and coefficient
			inline static polynomial<Type> monomial(Type c, unsigned int order) {
				
				polynomial m;
				m.coeff = std::vector<Type>(order + 1, Type(0));
				m.coeff[order] = c;
				return m;
			}


			// Friend operators to enable equations of the form
			// (Type) op. (polynomial<Type>)

			inline friend polynomial<Type> operator+(
				Type r, const polynomial<Type>& z) {

				return z + polynomial(r);
			}

			inline friend polynomial<Type> operator-(
				Type r, const polynomial<Type>& z) {

				return (z * -1) + polynomial(r);
			}

			inline friend polynomial<Type> operator*(
				Type r, const polynomial<Type>& z) {

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


			/// Convert the polynomial to string representation.
			inline operator std::string() {
				return to_string();
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
