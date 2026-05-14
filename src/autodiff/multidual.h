#ifndef THEORETICA_MULTIDUAL_H
#define THEORETICA_MULTIDUAL_H


///
/// @file multidual.h Multidual numbers
///

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../algebra/vec.h"
#include "../algebra/mat.h"


namespace theoretica {

	/// 
	/// @class multidual
	/// Multidual number algebra for functions
	/// of the form \f$f: \mathbb{R}^n \rightarrow \mathbb{R}\f$
	///
	template<unsigned int N = 0>
	class multidual {
		
		public:
			
			/// The real part of the multidimensional dual number
			real a;

			/// The dual part of the multidimensional dual number
			/// as a real vector.
			vec<real, N> v;

			// The template argument of the vector type used
			static constexpr size_t size_argument = N;


			/// Construct a multidual number
			/// as \f$(0 + \vec 0)\f$
			multidual() : a(0) {}


			/// Construct a multidual number from
			/// a real number and an N dimensional vector
			multidual(real r, vec<real, N> u) : a(r), v(u) {}


			/// Construct a multidual number from
			/// a real number
			multidual(real r) : a(r), v(vec<real, N>()) {}

			~multidual() = default;


			/// Initialize a multidual number from a real number
			inline multidual& operator=(real x) {
				a = x;
				v = vec<real, N>();
				return *this;
			}


			/// Get the real part
			inline real Re() const {
				return a;
			}


			/// Access the real part
			inline real& Re() {
				return a;
			}


			/// Get the multidual part
			inline vec<real, N> Dual() const {
				return v;
			}


			/// Access the multidual part
			inline vec<real, N>& Dual() {
				return v;
			}


			/// Get the i-th element of the multidual part,
			/// corresponding to the i-th independent variable
			/// in automatic differentiation.
			inline real Dual(unsigned int i) const {
				return v[i];
			}


			/// Access the i-th element of the multidual part,
			/// corresponding to the i-th independent variable
			/// in automatic differentiation.
			inline real& Dual(unsigned int i) {
				return v[i];
			}


			/// Get the multidual conjugate
			inline multidual conjugate() const {
				return multidual(a, -v);
			}


			/// Get the inverse of a multidual number
			inline multidual inverse() const {

				if(a == 0) {
					TH_MATH_ERROR("multidual::inverse", 0, MathError::DivByZero);
					return multidual(nan(), vec<real, N>(N, nan()));
				}

				return multidual(1.0 / a, v * (-1 / (a * a)));
			}


			/// Identity (for consistency)
			inline multidual operator+() const {
				return multidual(a, v);
			}


			/// Sum two multidual numbers
			inline multidual operator+(const multidual& other) const {
				return multidual(a + other.a, v + other.v);
			}


			/// Sum a real number to a multidual number
			inline multidual operator+(real r) const {
				return multidual(a + r, v);
			}


			/// Get the opposite of a multidual number
			inline multidual operator-() const {
				return multidual(-a, -v);
			}


			/// Subtract two multidual numbers
			inline multidual operator-(const multidual& other) const {
				return multidual(a - other.a, v - other.v);
			}


			/// Subtract a real number from a multidual number
			inline multidual operator-(real r) const {
				return multidual(a - r, v);
			}


			/// Multiply two multidual numbers
			inline multidual operator*(const multidual& other) const {

				// Equal sizes
				if (v.size() == other.v.size()) {
					return multidual(a * other.a, other.v * a + v * other.a);
				}

				// Trivial scalar cases
				if (v.size() == 0) {
					return multidual(a * other.a, other.v * a);
				}

				if (other.v.size() == 0) {
					return multidual(a * other.a, v * other.a);
				}

				// Treat missing entries as zero
				const unsigned int new_size = th::max(v.size(), other.v.size());
				vec<real, N> res_v;
				res_v.resize(new_size);

				for (unsigned int i = 0; i < res_v.size(); ++i) {
					real lhs = (i < other.v.size()) ? other.v[i] : 0.0;
					real rhs = (i < v.size()) ? v[i] : 0.0;
					res_v[i] = lhs * a + rhs * other.a;
				}

				return multidual(a * other.a, res_v);
			}


			/// Multiply a multidual number by a real number
			inline multidual operator*(real r) const {
				return multidual(a * r, v * r);
			}


			/// Dual division
			inline multidual operator/(const multidual& other) const {

				if(abs(a) < MACH_EPSILON) {
					TH_MATH_ERROR("multidual::operator/", 0, MathError::DivByZero);
					return multidual(nan(), vec<real, N>(N, nan()));
				}

				const real res_a = a / other.a;
				const real denom = other.a * other.a;

				// Equal sizes
				if (v.size() == other.v.size()) {
					return multidual(res_a, (v * other.a - other.v * a) / denom);
				}

				// Trivial scalar cases
				if (v.size() == 0) {
					return multidual(res_a, other.v * (-a / denom));
				}

				if (other.v.size() == 0) {
					return multidual(res_a, v / other.a);
				}

				// Treat missing entries as zero
				const unsigned int new_size = th::max(v.size(), other.v.size());
				vec<real, N> res_v;
				res_v.resize(new_size);

				for (unsigned int i = 0; i < res_v.size(); ++i) {

					real lhs = (i < v.size()) ? v[i] : 0.0;
					real rhs = (i < other.v.size()) ? other.v[i] : 0.0;

					res_v[i] = (lhs * other.a - rhs * a) / denom;
				}

				return multidual(res_a, res_v);
			}

			/// Divide a multidual number by a real number
			inline multidual operator/(real r) const {
				return multidual(a / r, v / r);
			}


			/// Sum a real number to this one
			inline multidual& operator+=(const multidual& other) {

				a += other.a;
				v += other.v;
				return *this;
			}


			/// Sum a real number to this multidual number
			inline multidual& operator+=(real r) {

				a += r;
				return *this;
			}

			/// Subtract a real number from this one
			inline multidual& operator-=(const multidual& other) {

				a -= other.a;
				v -= other.v;
				return *this;
			}

			/// Subtract a real number from this multidual number
			inline multidual& operator-=(real r) {

				a -= r;
				return *this;
			}

			/// Multiply this multidual number by another one
			inline multidual& operator*=(const multidual& other) {

				const real a_old = a;

				// Equal sizes
				if (v.size() == other.v.size()) {
					a = a_old * other.a;
					v = other.v * a_old + v * other.a;
					return *this;
				}

				// Trivial scalar cases
				if (v.size() == 0) {
					a = a_old * other.a;
					v = other.v * a_old;
					return *this;
				}

				if (other.v.size() == 0) {
					a = a_old * other.a;
					v = v * other.a;
					return *this;
				}

				// Treat missing entries as zero
				const unsigned int new_size = th::max(v.size(), other.v.size());
				vec<real, N> new_v;
				new_v.resize(new_size);

				for (unsigned int i = 0; i < new_v.size(); ++i) {

					real lhs = (i < other.v.size()) ? other.v[i] : 0.0;
					real rhs = (i < v.size()) ? v[i] : 0.0;

					new_v[i] = lhs * a_old + rhs * other.a;
				}

				a = a_old * other.a;
				v = new_v;

				return *this;
			}

			/// Multiply this multidual number by a real number
			inline multidual& operator*=(real r) {

				a *= r;
				v *= r;
				return *this;
			}

			/// Divide this multidual number by another one
			inline multidual& operator/=(const multidual& other) {

				const real a_old = a;
				const real denom = other.a * other.a;

				// Equal sizes
				if (v.size() == other.v.size()) {
					a = a_old / other.a;
					v = (v * other.a - other.v * a_old) / denom;
					return *this;
				}

				// Trivial scalar cases
				if (v.size() == 0) {
					a = a_old / other.a;
					v = other.v * (-a_old / denom);
					return *this;
				}

				if (other.v.size() == 0) {
					a = a_old / other.a;
					v = v / other.a;
					return *this;
				}

				// Treat missing entries as zero
				const unsigned int new_size = th::max(v.size(), other.v.size());
				vec<real, N> new_v;
				new_v.resize(new_size);

				for (unsigned int i = 0; i < new_v.size(); ++i) {

					real lhs = (i < v.size()) ? v[i] : 0.0;
					real rhs = (i < other.v.size()) ? other.v[i] : 0.0;

					new_v[i] = (lhs * other.a - rhs * a_old) / denom;
				}

				a = a_old / other.a;
				v = new_v;

				return *this;
			}

			/// Divide a multidual number by a real number
			inline multidual& operator/=(real r) {

				if(r == 0) {
					TH_MATH_ERROR("multidual::operator/=", 0, MathError::DivByZero);
					a = nan();
					v = vec<real, N>(N, nan());
					return *this;
				}

				a /= r;
				v /= r;

				return *this;
			}


			/// Check whether two multidual numbers have the same
			/// real and multidual parts
			inline bool operator==(const multidual& other) {
				return (a == other.a) && (v == other.v);
			}


			/// Construct an N-dimensional vector of multidual numbers
			/// to be passed as argument to a multidual function
			/// @param x A vector of real numbers containing the variables to pass
			template<typename Vector = vec<real, N>, enable_vector<Vector> = true>
			inline static vec<multidual<N>, N> make_argument(const Vector& x) {

				vec<multidual<N>, N> arg;
				arg.resize(x.size());

				for (unsigned int i = 0; i < x.size(); ++i)
					arg[i] = multidual<N>(
						x[i], vec<real, N>::euclidean_base(i, x.size())
					);

				return arg;
			}


			/// Extract the real vector from a vector of
			/// multidual numbers as a vec<real, N>
			inline static vec<real, N> extract_real(
				const vec<multidual<N>, N>& v) {

				vec<real, N> x;
				x.resize(v.size());

				for (unsigned int i = 0; i < N; ++i)
					x[i] = v[i].Re();

				return x;
			}


			/// Extract the dual matrix (Jacobian)
			/// from a vector of multidual numbers
			/// as a mat<N, N>
			inline static mat<real, N, N> extract_dual(
				const vec<multidual<N>, N>& v) {

				mat<real, N, N> J;
				J.resize(v.size(), v.size());

				for (unsigned int i = 0; i < N; ++i)
					for (unsigned int j = 0; j < N; ++j)
						J(j, i) = v[j].Dual(i);
				
				return J;
			}


			/// Extract the real vector and dual matrix
			/// from a vector of multidual numbers
			/// as a vec<real, N> and mat<real, N, N>
			inline static void extract(
				const vec<multidual<N>, N>& v,
				vec<real, N>& x,
				mat<real, N, N>& J) {

				for (unsigned int i = 0; i < N; ++i) {
					
					for (unsigned int j = 0; j < N; ++j)
						J(j, i) = v[j].Dual(i);
					
					x[i] = v[i].Re();
				}
			}


			/// Get the number of independent variables associated
			/// with the multidual number.
			unsigned int size() const {
				return v.size();
			}


			/// Change the size of the dual part of the number
			/// (only for dynamically allocated vectors)
			inline void resize(unsigned int size) {
				v.resize(size);
			}


			// Friend operators to enable equations of the form
			// (real) op. (multidual)
			
			inline friend multidual operator+(real a, const multidual& d) {
				return d + a;
			}

			inline friend multidual operator-(real a, const multidual& d) {
				return -d + a;
			}

			inline friend multidual operator*(real a, const multidual& d) {
				return d * a;
			}

			inline friend multidual operator/(real a, const multidual& d) {
				return a * d.inverse();
			}


#ifndef THEORETICA_NO_PRINT

			/// Convert the multidual number to string representation
			/// @param epsilon The character to use to represent epsilon
			inline std::string to_string(const std::string& epsilon = "ε") const {

				std::stringstream res;
				res << a << " + " << v << epsilon;

				return res.str();
			}


			/// Convert the multidual number to string representation.
			inline operator std::string() {
				return to_string();
			}


			/// Stream the multidual number in string representation
			/// to an output stream (std::ostream)
			inline friend std::ostream& operator<<(std::ostream& out, const multidual& obj) {
				return out << obj.to_string();
			}

#endif
		
	};

}

#endif
