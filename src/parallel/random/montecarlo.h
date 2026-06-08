
///
/// @file montecarlo.h Parallelized implementations of Monte Carlo methods.
///

#ifndef THEORETICA_PARALLEL_MONTECARLO_H
#define THEORETICA_PARALLEL_MONTECARLO_H

#include "algebra/algebra_types.h"
#include "core/core_traits.h"
#include "random/random.h"
#include "random/stoch_result.h"
#include "statistics/running.h"

#include <omp.h>
#include <functional>


namespace theoretica {
namespace random {
namespace parallel {


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

		// Statistics are computed online using Welford's method.
		stats::RunningMoments2<real> stats;

		// Pre-generate seeds for all threads to avoid contention
		int num_threads;
		#pragma omp parallel
		{
			#pragma omp single
			num_threads = omp_get_num_threads();
		}
		
		vec<uint64_t> seeds (num_threads);
		for (int i = 0; i < num_threads; ++i)
				seeds[i] = g();

		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			stats::RunningMoments2<real> local_stats;
			PRNG local_g = PRNG(seeds[id]);

			#pragma omp for nowait
			for (unsigned int i = 0; i < N; ++i) {
				const real x = random::uniform(a, b, local_g);
				local_stats.insert(f(x));
			}

			#pragma omp critical
			{
                // Combine statistics accross threads using Chan's formula
				stats.combine(local_stats);
			}
		}

		const real volume = (b - a);

		return stoch_result<real>{
			volume * stats.mean(),
			volume * sqrt(stats.variance() * (1.0 / stats.number())),
			stats.number()
		};
	}


    /// Approximate an integral by using Crude Monte Carlo integration.
	///
	/// @param f The multivariate function to integrate (\f$f: \mathbb{R}^n \rightarrow \mathbb{R}^m\f$)
	/// @param extremes A vector of the extremes of integration (e.g. a vector of vec2)
	/// @param g An already initialized PRNG
	/// @param N The number of points to sample (defaults to 1M)
	/// @tparam ThreadPRNG The type of the PRNG to use in each thread (defaults to PRNG).
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

		// Statistics are computed online using Welford's method.
		stats::RunningMoments2<ReturnType> stats;

		// Pre-generate seeds for all threads to avoid contention
		int num_threads;
		#pragma omp parallel
		{
			#pragma omp single
			num_threads = omp_get_num_threads();
		}
		
		vec<uint64_t> seeds (num_threads);
		for (int i = 0; i < num_threads; ++i)
				seeds[i] = g();

		#pragma omp parallel
		{
			int id = omp_get_thread_num();
			stats::RunningMoments2<ReturnType> local_stats;
			PRNG local_g = PRNG(seeds[id]);

			Vector v;
			v.resize(dim);

			#pragma omp for nowait
			for (unsigned int i = 0; i < N; ++i) {
				
				for (unsigned int k = 0; k < v.size(); ++k)
					v[k] = random::uniform(extremes[k][0], extremes[k][1], local_g);

				local_stats.insert(f(v));
			}

			#pragma omp critical
			{
                // Combine statistics accross threads using Chan's formula
				stats.combine(local_stats);
			}
		}

		real volume = 1;
		for (unsigned int i = 0; i < dim; ++i)
			volume *= (extremes[i][1] - extremes[i][0]);

		return stoch_result<ReturnType>{
			volume * stats.mean(),
			volume * sqrt(stats.variance() * (1.0 / stats.number())),
			stats.number()
		};
	}


    /// Approximate an integral by using the randomized Crude
	/// Quasi-Monte Carlo integration by sampling from the Weyl sequence.
	/// The total number of points N is distributed over N_shift random shifts,
	/// and the variance of the estimates is used to compute the final error.
	///
	/// @param f The function to integrate
	/// @param a The lower extreme of integration
	/// @param b The upper extreme of integration
	/// @param N The total number of points to generate (assumed to be divisible by N_shift)
	/// @param N_shift The number of random shifts to use
	template<typename RealFunction, typename PRNG>
	inline stoch_result<real> integral_qmc(
		RealFunction f, real a, real b,
		PRNG& g, size_t N = 1'000'000, size_t N_shift = 10
	) {

		stats::RunningMoments2<real> stats;

		// Pre-generate shifts to avoid contention during parallel region
		vec<real> shifts (N_shift);
		for (unsigned int i = 0; i < N_shift; ++i)
			shifts[i] = random::uniform(0, 1, g);

        #pragma omp parallel for
		for (unsigned int j = 0; j < N_shift; ++j) {

			// Random shift for the Weyl sequence
			real shift = shifts[j];

			// Distribute points over shifts
			const size_t size = N / N_shift;
			real sum_y = 0;

			for (unsigned int i = 0; i < size; ++i) {

				const real x = a + fract(i * INVPHI + shift) * (b - a);
				sum_y += f(x);
			}

            #pragma omp critical
            {	
				const real estimate = (b - a) * sum_y / size;
                stats.insert(estimate);
            }
		}
	
		return stoch_result<real>{
			stats.mean(),
			sqrt(stats.variance()),
			N
		};
	}

}}}

#endif
