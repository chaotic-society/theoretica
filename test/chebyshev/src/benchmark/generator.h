
///
/// @file generator.h Input generators for benchmarks.
///

#ifndef CHEYBYSHEV_GENERATOR_H
#define CHEYBYSHEV_GENERATOR_H

#include <functional>
#include "../core/common.h"
#include "../core/random.h"


namespace chebyshev {
namespace benchmark {

	namespace generator {

		/// Uniform generator over a domain
		inline auto uniform1D(long double a, long double b) {

			return [=](unsigned int i) {
				return random::uniform(a, b);
			};
		}

	}

}}

#endif
