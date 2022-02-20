#ifndef UROBORO_POLYNOMIAL_H
#define UROBORO_POLYNOMIAL_H

#include "./vec_buff.h"
#include "./real_analysis.h"


namespace uroboro {

	// A polynomial of arbitrary order with real coefficients
	template<typename T = real>
	class polynomial {
		public:
			std::vector<T> coeff;

			polynomial() : coeff() {}

			polynomial(const std::vector<T>& c) : coeff(c) {}
			
			~polynomial() {}

			polynomial(std::initializer_list<real> l) : coeff(l) {}


			// Access i-th coefficient
			inline T& operator[](int i) {
				return coeff[i];
			}


			// Get the i-th by value
			inline T get(int i) const {
				return coeff[i];
			}


			// Evaluate the polynomial using x as variable
			inline T eval(T x) const {

				T res = 0;

				// TO-DO Implement Horner's method

				for (int i = 0; i < coeff.size(); ++i) {
					if(coeff[i] != 0) {
						res += coeff[i] * uroboro::pow(x, i);
					}
				}

				return res;
			}


			// Evaluate the polynomial using x as variable
			inline T operator()(T x) const {
				return eval(x);
			}


			// Return the nth order coefficient
			inline T& operator[](unsigned int i) {
				return coeff[i];
			}


			// Find the true order of the polynomial (ignoring null coefficients)
			inline int find_order() const {

				for (int i = coeff.size() - 1; i >= 0; --i) {
					if(coeff[i] != 0)
						return i;
				}

				return 0;
			}


			// Remove higher order null coefficients
			inline void trim() {

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


			// TO-DO Polynomial division


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


			// Check whether two polynomials are equal
			inline bool operator==(const polynomial& other) const {

				for (int i = 0; i < min(other.size(), this->size()); ++i) {
					if(other.coeff[i] != coeff[i])
						return false;
				}

				if(size() > other.size()) {
					for (int i = 0; i < size(); ++i) {
						if(coeff[i] != 0)
							return false;
					}
				} else {
					for (int i = 0; i < other.size(); ++i) {
						if(other.coeff[i] != 0)
							return false;
					}
				}

				return true;
			}
		
	};
	
}

#endif
