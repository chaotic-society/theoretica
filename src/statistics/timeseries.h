#ifndef THEORETICA_TIMESERIES_H
#define THEORETICA_TIMESERIES_H

#include "../core/constants.h"


namespace theoretica {

	namespace timeseries {


		// Compute the lag-n autocorrelation of a times series
		template<typename Dataset>
		inline real autocorrelation(const Dataset& x, unsigned int lag = 1) {

			if(lag > x.size()) {
				TH_MATH_ERROR("th::timeseries::autocorrelation",
					lag, INVALID_ARGUMENT);
				return nan();
			}

			real res = 0;

			for (int i = 0; i < x.size() - lag; ++i)
				res += x[i] * x[i + lag];

			return res;
		}


	}
}

#endif
