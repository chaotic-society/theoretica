
///
/// @file generator.h Input generators for benchmarks.
///

#ifndef CHEYBYSHEV_GENERATOR_H
#define CHEYBYSHEV_GENERATOR_H

#include <functional>
#include "../core/common.h"


namespace chebyshev {
namespace benchmark {

	namespace generator {

		/// Uniform generator over a domain
		template<typename FloatType = long double>
		inline auto uniform1D(long double a, long double b) {

			return [=](unsigned int i) {
				return random_uniform(a, b);
			};
		}

	}

}}

#endif
