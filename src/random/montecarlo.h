
///
/// @file montecarlo.h Monte Carlo methods.
///

#ifndef THEORETICA_MONTECARLO_H
#define THEORETICA_MONTECARLO_H

#include "../core/core_traits.h"
#include "./random.h"
#include "./quasirandom.h"
#include "./sampling.h"
#include "./stoch_result.h"
#include "../statistics/running.h"
#include "../algebra/vec_functions.h"

#ifdef _OPENMP
#include "../parallel/random/montecarlo.h"
#endif


namespace theoretica {
namespace random {


	/// Approximate an integral by using Crude Monte Carlo integration.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param g An already initialized PRNG
	/// @param N The number of points to sample (defaults to 1M)
	template<typename RealFunction, typename PRNG>
	inline stoch_result<real> integral_mc(
		RealFunction f, real a, real b,
		PRNG& g, size_t N = 1'000'000
	) {

		// Parallel dispatch
		#ifdef _OPENMP
		if (N >= 10'000)
			return parallel::integral_mc(f, a, b, g, N);
		#endif

		// Statistics are computed online using Welford's method.
		stats::RunningMoments2<real> stats;
		for (unsigned int i = 0; i < N; ++i) {
			const real x = random::uniform(a, b, g);
			stats.insert(f(x));
		}

		return stoch_result<real>{
			stats.mean() * (b - a),
			sqrt(stats.variance() * ((b - a) * (b - a) / stats.number())),
			stats.number()
		};
	}


	/// Approximate a multivariate integral by using Crude Monte Carlo integration.
	///
	/// @param f The multivariate function to integrate (\f$f: \mathbb{R}^n \rightarrow \mathbb{R}^m\f$)
	/// @param extremes A vector of the extremes of integration (e.g. a vector of vec2)
	/// @param g An already initialized PRNG
	/// @param N The number of points to sample (defaults to 1M)
	template <
		typename Function,
		typename DomainVector = vec<vec2>,
		typename PRNG,
		typename ReturnType = return_type_t<Function>
	>
	inline stoch_result<ReturnType> integral_mc(
		Function f, DomainVector extremes,
		PRNG& g, size_t N = 1'000'000
	) {

		using Vector = typename _internal::func_helper<Function>::first_arg_type;
		const size_t dim = extremes.size();


		// Parallel dispatch
		#ifdef _OPENMP
		if (N * dim >= 10'000)
			return parallel::integral_mc(f, extremes, g, N);
		#endif


		// Statistics are computed online using Welford's method.
		stats::RunningMoments2<ReturnType> stats;

		// Initialize sum by evaluating a first point
		Vector v;
		v.resize(dim);

		// Sample the function at random points in the integration region
		for (unsigned int i = 0; i < N; ++i) {
			
			for (unsigned int k = 0; k < v.size(); ++k)
				v[k] = random::uniform(extremes[k][0], extremes[k][1], g);

			stats.insert(f(v));
		}

		// Volume of the integral domain
		real volume = 1;
		for (unsigned int i = 0; i < dim; ++i)
			volume *= (extremes[i][1] - extremes[i][0]);
		

		// Return result as a stochastic result
		return stoch_result<ReturnType>{
			stats.mean() * volume,
			volume * sqrt(stats.variance() * (1.0 / stats.number())),
			stats.number()
		};
	}


	/// Approximate an integral by using the randomized Crude
	/// Quasi-Monte Carlo integration by sampling from the Weyl sequence.
	/// The total number of points N is distributed over M random shifts,
	/// and the variance of the estimates is used to compute the final error.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param N The total number of points to generate (expected to be a multiple of M)
	/// @param M The number of random shifts to use
	template<typename RealFunction, typename PRNG>
	inline stoch_result<real> integral_qmc(
		RealFunction f, real a, real b,
		PRNG& g, size_t N = 1'000'000, size_t M = 10
	) {

		// Parallel dispatch
		#ifdef _OPENMP
		if (N >= 1'000'000)
			return parallel::integral_qmc(f, a, b, g, N, M);
		#endif


		stats::RunningMoments2<real> stats;
		for (unsigned int run = 0; run < M; ++run) {

			// Random shift for the Weyl sequence
			const real shift = random::uniform(0, 1, g);

			// Distribute points over shifts
			const size_t size = N / M;

			real sum_y = 0;
			for (unsigned int i = 0; i < size; ++i) {

				const real x = a + fract(i * INVPHI + shift) * (b - a);
				sum_y += f(x);
			}

			const real estimate = (b - a) * sum_y / size;
			stats.insert(estimate);
		}
	
		return stoch_result<real>{
			stats.mean(),
			sqrt(stats.variance()),
			N
		};
	}


	/// Approximate a multivariate integral by using the randomized Crude
	/// Quasi-Monte Carlo integration by sampling from the Weyl sequence.
	/// The total number of points N is distributed over M random shifts,
	/// and the variance of the estimates is used to compute the final error.
	///
	/// @param f The function to integrate
	/// @param extremes A vector of the extremes of integration
	/// @param alpha An irrational number
	/// @param N The number of points to generate (expected to be a multiple of M)
	/// @param M The number of random shifts to use
	template <
		typename Function,
		typename DomainVector = vec<vec2>,
		typename PRNG,
		typename ReturnType = return_type_t<Function>,
		typename AlphaVector = vec<real>
	>
	inline stoch_result<ReturnType> integral_qmc(
		Function f, DomainVector extremes,
		PRNG& g, size_t N = 1'000'000, size_t M = 10,
		AlphaVector alphas = {}
	) {

		using Vector = typename _internal::func_helper<Function>::first_arg_type;
		const size_t dim = extremes.size();

		// If not provided, compute the optimal
		// parameters for the Weyl sequence.
		if (alphas.size() == 0) {

			alphas.resize(dim);

			// Find the only positive root of the polynomial: x^s+1 - x - 1 = 0
			const real base = 1.0 / root_bisect([=](real x) {
				return pow(x, dim + 1) - x - 1;
			}, 0.0, 2.0);

			// Optimal parameters are powers of the base
			real alpha = base;
			for (unsigned int i = 0; i < dim; ++i) {
				alphas[i] = alpha;
				alpha *= base;
			}
		}

		// Volume of the integral domain
		real volume = 1.0;
		for (unsigned int i = 0; i < dim; ++i)
			volume *= (extremes[i][1] - extremes[i][0]);


		stats::RunningMoments2<ReturnType> stats;
		for (unsigned int run = 0; run < M; ++run) {

			// Random shift for the Weyl sequence
			Vector shift;
			shift.resize(dim);
			for (unsigned int i = 0; i < shift.size(); ++i)
				shift[i] = random::uniform(0, 1, g);

			// Distribute points over shifts
			const size_t size = N / M;

			Vector x;
			x.resize(dim);

			ReturnType sum_y = 0;
			for (unsigned int i = 0; i < size; ++i) {

				for (unsigned int k = 0; k < x.size(); ++k) {
					const real xi = fract(i * alphas[k] + shift[k]);
					x[k] = extremes[k][0] + xi * (extremes[k][1] - extremes[k][0]);
				}

				sum_y += f(x);
			}

			stats.insert(volume * sum_y / size);
		}

		return stoch_result<ReturnType>{
			stats.mean(),
			sqrt(stats.variance()),
			N
		};
	}


	/// Approximate an integral using Crude Monte Carlo integration with
	/// importance sampling.
	///
	/// @param f The function to integrate
	/// @param g The importance function (normalized)
	/// @param g_sampler A functor which samples the importance
	/// distribution (e.g. a PdfSampler<PRNG> object).
	/// @param N The number of points to generate
	template <
		typename RealFunction,
		typename PdfFunction,
		typename PdfSampler
	>
	inline real integral_impsamp(
		RealFunction f, PdfFunction g,
		PdfSampler g_sampler, size_t N = 1'000'000) {

		stats::RunningMoments2<real> stats;
		for (unsigned int i = 0; i < N; ++i) {
			const auto v = g_sampler();
			stats.insert(f(v) / g(v));
		}		

		return stoch_result<real>{
			stats.mean(),
			sqrt(stats.variance()),
			stats.number()
		};
	}


	/// Generate a Monte Carlo sample of values of a given real function
	/// of arbitrary variables following the given distributions.
	///
	/// @tparam SampleVector The type of the vector to return (defaults to vec<real>)
	/// @param N The size of the sample.
	/// @param f The function to sample, taking as arguments the random variables,
	/// in the same order as the samplers.
	/// @param samplers A variadic list of random variable samplers, which are
	/// functors that return a random value following the desired distribution
	/// (e.g. a PdfSampler<PRNG> object).
	/// @return The sampled function values.
	template <
		typename SampleVector = vec<real>,
		typename Function,
		typename ...PdfSamplers
	>
	SampleVector sample_function(size_t N, Function f, PdfSamplers... samplers) {

		SampleVector sample;
		sample.resize(N);

		for (unsigned int i = 0; i < N; ++i)
			sample[i] = f(samplers()...);

		return sample;
	}


	/// Generate a Monte Carlo sample of values of a given real function
	/// of arbitrary variables following the given distributions,
	/// overwriting every element of the given vector.
	///
	/// @tparam SampleVector The type of the vector to return (defaults to vec<real>)
	/// @param v The vector to overwrite with sampled values.
	/// @param f The function to sample, taking as arguments the random variables,
	/// in the same order as the samplers.
	/// @param samplers A variadic list of random variable samplers, which are
	/// functors that return a random value following the desired distribution
	/// (e.g. a PdfSampler<PRNG> object).
	/// @return The sampled function values.
	template <
		typename Vector,
		typename Function,
		typename ...PdfSamplers
	>
	Vector& sample_function(Vector& v, Function f, PdfSamplers... samplers) {

		for (unsigned int i = 0; i < v.size(); ++i)
			v[i] = f(samplers()...);

		return v;
	}

}}


#endif
