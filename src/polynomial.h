#ifndef UROBORO_POLYNOMIAL_H
#define UROBORO_POLYNOMIAL_H

#include "./vec_buff.h"
#include "./common.h"


namespace uroboro {

	class polynomial {
		public:
			vec_buff coeff;

			polynomial() : coeff() {}

			polynomial(const vec_buff& c) : coeff(c) {}
			
			~polynomial() {}

			polynomial(std::initializer_list<real> l) : coeff(l) {}


			// Access i-th coefficient
			inline real& operator[](int i) {
				return coeff[i];
			}


			inline real get(int i) const {
				return coeff[i];
			}


			// Evaluate the polynomial using x as variable
			inline real eval(real x) const {

				real res = 0;

				// TO-DO Implement Horner's method

				for (int i = 0; i < coeff.size(); ++i) {
					if(coeff[i] != 0) {
						res += coeff[i] * uroboro::pow(x, i);
					}
				}

				return res;
			}


			// Evaluate the polynomial using x as variable
			inline real operator()(real x) const {
				return eval(x);
			}


			// Find the true order of the polynomial (without counting null coefficients)
			inline int find_order() const {

				for (int i = coeff.size() - 1; i >= 0; --i) {
					if(coeff[i] != 0)
						return i;
				}

				return 0;
			}


			// Remove higher order null coefficients
			inline int trim() {

				for (int i = coeff.size() - 1; i >= 0; --i) {
					if(coeff[i] != 0)
						break;

					coeff.pop_back();
				}
			}


			// Get the number of coefficients
			inline size_t size() const {
				return coeff.size();
			}


			// Sum two polynomials
			inline polynomial operator+(const polynomial& p) const {

				polynomial r(coeff);

				for (int i = 0; i < uroboro::min(r.size(), p.size()); ++i) {
					r[i] += p.get(i);
				}

				return r;
			}


			// Subtract a polynomial from another
			inline polynomial operator-(const polynomial& p) const {
				
				polynomial r(coeff);

				for (int i = 0; i < uroboro::min(r.size(), p.size()); ++i) {
					r[i] -= p.get(i);
				}

				return r;
			}


			// Multiply two polynomials
			inline polynomial operator*(const polynomial& p) const {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1);

				for (int i = 0; i < size(); ++i) {
					for (int j = 0; j < p.size(); ++j) {
						r[i + j] += coeff[i] * p.get(j);
					}
				}

				return r;
			}


			// Sum a polynomial to this one
			inline polynomial& operator+=(const polynomial& p) {

				for (int i = 0; i < uroboro::min(size(), p.size()); ++i) {
					coeff[i] += p.get(i);
				}

				return *this;
			}


			// Subtract a polynomial from this one
			inline polynomial& operator-=(const polynomial& p) {

				for (int i = 0; i < uroboro::min(size(), p.size()); ++i) {
					coeff[i] -= p.get(i);
				}

				return *this;
			}


			// Multiply two polynomials
			inline polynomial& operator*=(const polynomial& p) {

				polynomial r = polynomial();
				r.coeff.resize(this->size() + p.size() - 1);

				for (int i = 0; i < size(); ++i) {
					for (int j = 0; j < p.size(); ++j) {
						r[i + j] += coeff[i] * p.get(j);
					}
				}

				*this = r;

				return *this;
			}
		
	};
	
}

#endif
