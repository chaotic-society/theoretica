#ifndef THEORETICA_REGRESSION_H
#define THEORETICA_REGRESSION_H

namespace theoretica {

	/// Fit a set of data points to a model.
	///
	/// @param data The set of (x, y) data points
	/// @param f_model The model of the regression, taking in
	/// an x value and a list of the parameters.
	/// @param guess The initial guess for the parameters,
	/// defaults to {1, ..., 1}.
	template<unsigned int N>
	inline vec<N> fit(
		std::vector<vec<2>> data,
		multidual<N>(*f_model)(real, vec<N, multidual<N>>),
		vec<N> guess = vec<N>(1)) {

		return minimize<N>(
			[f_model, data](vec<N, multidual<N>> param) {

				multidual<N> loss = multidual<N>(0);

				for (unsigned int i = 0; i < data.size(); ++i)
					loss += square(
						data.at(i).get(1) - f_model(data.at(i).get(0), param)
					);

				return loss / data.size();
			}, guess);
	}


}


#endif
