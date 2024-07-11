
///
/// @file estimator.h Default precision estimators.
///

#ifndef CHEBYSHEV_ESTIMATOR
#define CHEBYSHEV_ESTIMATOR


#include <functional>
#include <cmath>

#include "../core/common.h"
#include "../core/random.h"
#include "./prec_structures.h"


namespace chebyshev {
namespace prec {


	/// @namespace chebyshev::prec::estimator Precision estimators.
	namespace estimator {


		/// Use Simpson's quadrature scheme to approximate error integrals
		/// for univariate real functions (endofunctions on real number types).
		/// The estimator is returned as a lambda function.
		template<typename FloatType = double>
		inline auto quadrature1D() {

			return [](
				EndoFunction<FloatType> funcApprox,
				EndoFunction<FloatType> funcExpected,
				estimate_options<FloatType, FloatType> options) {

				if(options.domain.size() != 1)
					throw std::runtime_error(
						"estimator::quadrature1D only works on mono-dimensional domains");

				interval domain = options.domain[0];

				FloatType sum = 0;
				FloatType sum_sqr = 0;
				FloatType sum_abs = 0;
				FloatType max = 0;

				const FloatType length = domain.length();
				const FloatType dx = length / options.iterations;
				FloatType x;
				FloatType coeff;

				FloatType diff = std::abs(funcApprox(domain.a) - funcExpected(domain.a));

				sum += diff;
				sum_sqr += diff * diff;
				sum_abs += std::abs(funcExpected(domain.a));
				max = diff;

				for (unsigned int i = 1; i < options.iterations; ++i) {

					x = domain.a + i * dx;
					diff = std::abs(funcApprox(x) - funcExpected(x));

					if(diff > max)
						max = diff;

					if(i % 2 == 0)
						coeff = 2;
					else
						coeff = 4;

					sum += coeff * diff;
					sum_sqr += coeff * diff * diff;
					sum_abs += coeff * funcExpected(x);
				}

				diff = std::abs(funcApprox(domain.b) - funcExpected(domain.b));

				sum += diff;
				sum_sqr += diff * diff;
				sum_abs += std::abs(funcExpected(domain.b));
				
				if(diff > max)
					max = diff;

				estimate_result res {};
				res.absErr = sum;
				res.maxErr = max;
				res.meanErr = (sum * dx / 3.0) / length;
				res.rmsErr = std::sqrt((sum_sqr * dx / 3.0) / length);
				res.relErr = std::abs((sum * dx / 3.0) / (sum_abs * dx / 3.0));
				
				return res;
			};
		}


		/// Use crude Monte Carlo integration to approximate error integrals
		/// for univariate real functions.
		template<typename FloatType = double>
		inline auto montecarlo1D() {

			// Return an 1-dimensional Monte Carlo estimator
			// as a lambda function
			return [](
				EndoFunction<FloatType> funcApprox,
				EndoFunction<FloatType> fExpected,
				estimate_options<FloatType, FloatType> options) {

				if(options.domain.size() != 1)
					throw std::runtime_error(
						"The estimation domain's dimension must be 1"
						"in prec::estimator::montecarlo1D");

				FloatType sum = 0;
				FloatType sum_sqr = 0;
				FloatType sum_abs = 0;
				FloatType max = 0;
				const FloatType length = options.domain[0].length();

				for (int i = 0; i < options.iterations; ++i) {
					
					FloatType x = random::uniform(options.domain[0].a, options.domain[0].b);

					const FloatType diff = std::abs(funcApprox(x) - funcExpected(x));

					if(max < diff)
						max = diff;

					sum += diff;
					sum_sqr += diff * diff;
					sum_abs += std::abs(funcExpected(x));
				}

				estimate_result res {};
				res.maxErr = max;
				res.meanErr = sum / options.iterations;
				res.absErr = sum * (length / options.iterations);
				res.rmsErr = sum_sqr * (length / options.iterations);
				res.relErr = sum / sum_abs;

				return res;
			};
		}


		/// Use crude Monte Carlo integration to approximate error integrals
		/// for multivariate real functions.
		///
		/// @param dimensions The dimension of the space of inputs
		/// @note You may specify a custom vector type to use as input,
		/// but it must provide a constructor taking in the number of elements.
		template<typename FloatType = double, typename Vector = std::vector<FloatType>>
		inline auto montecarlo(unsigned int dimensions) {

			// Return an n-dimensional Monte Carlo estimator
			// as a lambda function
			return [dimensions](
				std::function<FloatType(Vector)> funcApprox,
				std::function<FloatType(Vector)> fExpected,
				estimate_options<FloatType, FloatType> options) {

				if(options.domain.size() != dimensions)
					throw std::runtime_error(
						"The estimation domain's dimension does not match "
						"the instantiated number of dimensions "
						"in estimator::montecarlo");

				FloatType sum = 0;
				FloatType sum_sqr = 0;
				FloatType sum_abs = 0;
				FloatType max = 0;

				// Compute the measure of a multi-interval
				FloatType volume = 1;
				for (interval k : options.domain)
					volume *= k.length();

				Vector x (dimensions);

				for (int i = 0; i < options.iterations; ++i) {
					
					random::sample_uniform(x, options.domain);

					const FloatType diff = std::abs(funcApprox(x) - funcExpected(x));

					if(max < diff)
						max = diff;

					sum += diff;
					sum_sqr += diff * diff;
					sum_abs += std::abs(funcExpected(x));
				}

				estimate_result res {};
				res.maxErr = max;
				res.meanErr = sum / options.iterations;
				res.absErr = sum * (volume / options.iterations);
				res.rmsErr = sum_sqr * (volume / options.iterations);
				res.relErr = sum / sum_abs;

				return res;
			};
		}


	}

}}


#endif
