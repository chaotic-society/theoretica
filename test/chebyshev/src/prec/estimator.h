
///
/// @file estimator.h Default precision estimators.
///

#ifndef CHEBYSHEV_ESTIMATOR
#define CHEBYSHEV_ESTIMATOR


#include <functional>
#include <cmath>

#include "../core/common.h"
#include "./prec_structures.h"


namespace chebyshev {
namespace prec {


	/// @namespace chebyshev::prec::estimator Precision estimators.
	namespace estimator {


		/// Use Simpson's quadrature scheme to approximate error integrals
		/// for univariate real functions.
		template<typename FloatType = double>
		estimate_result quadrature1D(
			RealFunction<FloatType> funcApprox,
			RealFunction<FloatType> funcExpected,
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
		}


		/// Use crude Monte Carlo integration to approximate error integrals
		/// for multivariate real functions.
		template<typename FloatType = double>
		estimate_result montecarlo(unsigned int dimensions) {

			// Return an n-dimensional Monte Carlo estimator
			// as a lambda function
			return [dimensions](
				std::function<FloatType(std::vector<FloatType>)> funcApprox,
				std::function<FloatType(std::vector<FloatType>)> fExpected,
				estimate_options<FloatType, FloatType> options) {

				if(options.domain.size() != dimensions)
					throw std::runtime_error(
						"The estimation domain's dimension does not match "
						"the instantiated number of dimensions "
						"in estimator::montecarlo");

				FloatType sum = 0;
				FloatType sum_sqr = 0;
				FloatType sum_abs = 0;
				FloatType max = -std::numeric_limits<FloatType>::infinity();

				// Compute the measure of a multi-interval
				FloatType volume = 1;
				for (interval k : options.domain)
					volume *= k.length();

				std::vector<FloatType> x (dimensions);

				for (int i = 0; i < options.iterations; ++i) {
					
					sample_uniform(x, options.domain);

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
