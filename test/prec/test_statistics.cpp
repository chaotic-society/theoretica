
/// @file test_statistics.cpp Test cases for statistical functions.

#include "theoretica.h"
#include "chebyshev.h"

using namespace chebyshev;
using namespace theoretica;


// Sample statistics stochastic estimator.
//
// Given a sample v with a certain known statistic f(v), the estimator
// generates N derived samples from v which are shifted by a random quantity each.
// The test case requires passing a function which computes the expected statistic
// of the sample, given the random shift. This process is repeated over N samples,
// deriving an estimate of the numerical error.
template<typename Sample, typename Statistic, typename ShiftStatistic>
auto estimate_stat(
	prec::prec_context& ctx, const std::string& name,
	Sample v, Statistic stat, ShiftStatistic shifted) {

	using Result = typename std::decay<decltype(stat(v))>::type;
	using Options = prec::estimate_options<Result, Sample>;

	// Construct estimator in-place.
	auto estimator = [v, stat, shifted, &ctx, name](
		std::function<Result(Sample)> statFunc,
		std::function<Result(Sample)> expectedFunc,
		Options opt) {

		auto rnd = ctx.random->get_rnd();
		long double sum = 0.0;
		long double sum2 = 0.0;
		long double max = -std::numeric_limits<long double>::infinity();

		for (unsigned int i = 0; i < opt.iterations; ++i) {

			auto sample = v;
			real shift = rnd.uniform(opt.domain[0].a, opt.domain[0].b);
			for (auto& x : sample)
				x += shift;
			
			auto computed = statFunc(sample);
			auto expected = shifted(shift);
			(void) expectedFunc;

			real diff = std::abs(computed - expected);
			sum += diff;
			sum2 += diff * diff;

			if (diff > max)
				max = diff;
		}

		prec::estimate_result res;
			res.name = name;
		res.absErr = ch::get_nan<long double>();
		res.meanErr = sum / opt.iterations;
		res.rmsErr = std::sqrt(sum2 / opt.iterations);
		res.maxErr = max;
		res.tolerance = opt.tolerance;
		res.failed = opt.fail(res);
		return res;
	};

	auto opt = Options(
		prec::interval(-1E+06, 1E+06), estimator
	);

	ctx.estimate(name, stat, [shifted](vec<real>) { return shifted(0.0); }, opt);
}


int main(int argc, char const *argv[]) {

	auto ctx = prec::make_context("statistics", argc, argv);
	ctx.output->settings.outputFiles = {"test/prec/prec_statistics.csv"};
	ctx.settings.defaultIterations = 10'000;
	ctx.settings.defaultTolerance = 1E-08;
	
	prec::equation_options<vec<real>> vec_opt (1E-08, prec::distance::euclidean<vec<real>>);

	const vec<real> data = {2, 3, 5, 4, 7, 6, 10, 8, 12, 9};
	const vec<real> data2 = {4, 5, 11, 8, 14, 13, 21, 17, 24, 19};
	const vec<real> hist_data = {28, 31, 35, 29, 42, 37, 33, 48, 30, 39, 45, 34};
	const vec<real> autocorrelation_data = {1, 2, 3, 4};

	const real mean1 = 66 / 10.0;
	const real variance1 = 154.0 / 15.0;
	const real covariance1 = 967.0 / 45.0;

	// statistics.h
	{

		// Test statistical functions with stochastic shifts on the sample.
		estimate_stat(
			ctx, "stats::mean", data,
			stats::mean<vec<real>>,
			[&mean1](real shift) { return mean1 + shift; }
		);

		estimate_stat(
			ctx, "stats::variance", data,
			[](const vec<real>& v) { return stats::variance(v); },
			[&variance1](real shift) { return variance1; }
		);

		estimate_stat(
			ctx, "stats::stdev", data,
			[](const vec<real>& v) { return stats::stdev(v); },
			[&variance1](real shift) { return std::sqrt(variance1); }
		);


		// Evaluate statistics for a known, trivial sample, preventing regressions.

		ctx.equals("stats::mean", stats::mean(data), mean1);
		ctx.equals("stats::range", stats::range(data), 10.0);
		ctx.equals("stats::semidispersion", stats::semidispersion(data), 5.0);
		ctx.equals("stats::total_sum_squares", stats::total_sum_squares(data), 92.4);
		ctx.equals("stats::variance", stats::variance(data), variance1);
		ctx.equals("stats::variance (population)", stats::variance(data, 0), 9.24);
		ctx.equals("stats::stdev", stats::stdev(data), th::sqrt(variance1));
		ctx.equals("stats::stdom", stats::stdom(data), th::sqrt(154.0 / 15.0 / data.size()));

		{
			real moments_mean = 0;
			real moments_variance = 0;
			stats::moments2(data, moments_mean, moments_variance);
			ctx.equals("stats::moments2 (mean)", moments_mean, mean1);
			ctx.equals("stats::moments2 (variance)", moments_variance, variance1);
		}

		ctx.equals("stats::standard_relative_error",
			stats::standard_relative_error(data), th::sqrt(variance1) / mean1
		);

		ctx.equals("stats::covariance",
			stats::covariance(data, data2), covariance1
		);

		ctx.equals("stats::correlation_coefficient",
			stats::correlation_coefficient(data, data2), 0.9955833328
		);
		
		ctx.equals("stats::autocorrelation", stats::autocorrelation(data), 0.5242424242);
		ctx.equals("stats::autocorrelation (lag 0)", stats::autocorrelation(autocorrelation_data, 0), 1.0);
		ctx.equals("stats::autocorrelation (lag 2)", stats::autocorrelation(autocorrelation_data, 2), -0.3);
		ctx.equals("stats::autocorrelation (invalid lag)",
			std::isnan(stats::autocorrelation(autocorrelation_data, autocorrelation_data.size())), true
		);
		ctx.equals("stats::absolute_deviation", stats::absolute_deviation(data), 2.6);
		ctx.equals("stats::skewness", stats::skewness(data), 0.1444548871);
		ctx.equals("stats::kurtosis", stats::kurtosis(data), -1.3975324675);

		ctx.equals(
			"stats::propagate_sum",
			stats::propagate_sum(vec<real>{3, 4}), 5.0
		);
		
		ctx.equals(
			"stats::propagate_product",
			stats::propagate_product(vec<real>{1, 2}, vec<real>{2, 4}),
			th::sqrt(0.5)
		);

		ctx.equals(
			"stats::gaussian_expectation",
			stats::gaussian_expectation([](real x) {return x * x * x;}, 0, 100),
			0.0
		);

		ctx.equals(
			"stats::gaussian_expectation",
			stats::gaussian_expectation([](real x) {return x * x;}, 2, 3),
			13.0
		);

		ctx.equals("stats::z_score", stats::z_score(8, 5, 1.5), 2.0);

		ctx.equals(
			"stats::normalize_z_score",
			stats::normalize_z_score(data),
			vec<real>{-4.6, -3.6, -1.6, -2.6, 0.4, -0.6, 3.4, 1.4, 5.4, 2.4} / th::sqrt(variance1),
			vec_opt
		);

		ctx.equals(
			"stats::chi_square",
			stats::chi_square(vec<real>{10, 20}, vec<real>{8, 22}, vec<real>{1, 2}), 5.0
		);

		ctx.equals(
			"stats::chi_square_linear",
			stats::chi_square_linear(data, data2, vec<real>(10, 1.0), 0, 2),
			6.0
		);

		ctx.equals("stats::total_sum_squares (empty)",
			std::isnan(stats::total_sum_squares(vec<real>{})), true
		);
		ctx.equals("stats::variance (singleton)",
			std::isnan(stats::variance(vec<real>{3})), true
		);
		ctx.equals("stats::covariance (mismatched)",
			std::isnan(stats::covariance(data, vec<real>{1})), true
		);
		ctx.equals("stats::propagate_product (mismatched)",
			std::isnan(stats::propagate_product(vec<real>{1}, vec<real>{1, 2})), true
		);
		ctx.equals("stats::propagate_product (zero mean)",
			std::isnan(stats::propagate_product(vec<real>{1}, vec<real>{0})), true
		);
		ctx.equals("stats::chi_square (mismatched)",
			std::isnan(stats::chi_square(vec<real>{1}, vec<real>{1, 2}, vec<real>{1})), true
		);
		ctx.equals("stats::chi_square (zero uncertainty)",
			std::isnan(stats::chi_square(vec<real>{1}, vec<real>{1}, vec<real>{0})), true
		);


		const stat_function gaussian_stat = [](real x, const vec<real>& theta) {
			return distribution::gaussian(x, theta);
		};
		
		const vec<real> gaussian_parameters = {0, 1};

		const real gaussian_likelihood =
			distribution::gaussian(2, 0, 1) * distribution::gaussian(3, 0, 1)
			* distribution::gaussian(5, 0, 1) * distribution::gaussian(4, 0, 1)
			* distribution::gaussian(7, 0, 1) * distribution::gaussian(6, 0, 1)
			* distribution::gaussian(10, 0, 1) * distribution::gaussian(8, 0, 1)
			* distribution::gaussian(12, 0, 1) * distribution::gaussian(9, 0, 1);

		ctx.equals("stats::likelihood",
			stats::likelihood(data, gaussian_parameters, gaussian_stat),
			gaussian_likelihood
		);

		ctx.equals("stats::log_likelihood",
			stats::log_likelihood(data, gaussian_parameters, gaussian_stat),
			th::ln(gaussian_likelihood)
		);


		// Check specific values of chi-squared p-value
		ctx.equals("stats::pvalue_chi2 (invalid ndf)", std::isnan(stats::pvalue_chi2(0, 0)), true);

		ctx.equals("stats::pvalue_chi2(0, 1)", stats::pvalue_chi2(0, 1), 1.0, 1E-8);
		ctx.equals("stats::pvalue_chi2(0, 2)", stats::pvalue_chi2(0, 2), 1.0, 1E-8);
		ctx.equals("stats::pvalue_chi2(2, 2)", stats::pvalue_chi2(2, 2), th::exp(-1.0), 1E-8);
		ctx.equals("stats::pvalue_chi2(3.8414588207, 1)", stats::pvalue_chi2(3.8414588207, 1), 0.05, 1E-8);
		ctx.equals("stats::pvalue_chi2(5.9914645471, 2)",stats::pvalue_chi2(5.9914645471, 2), 0.05, 1E-8);
	}


	// distributions.h
	{
		const auto real_opt = prec::estimate_options<real, real>(
			prec::interval(-10, 10), prec::estimator::quadrature1D()
		);

		const auto exp_opt = prec::estimate_options<real, real>(
			prec::interval(0.001, 15), prec::estimator::quadrature1D()
		);

		const auto beta_opt = prec::estimate_options<real, real>(
			prec::interval(0.001, 0.999), prec::estimator::quadrature1D()
		);


		// Check the functional form for constant parameters

		ctx.estimate("distribution::gaussian",
			[](real x) {return distribution::gaussian(x, 0, 1);},
			[](real x) {return th::exp(-x * x / 2) / SQRT2 / SQRTPI;},
			real_opt
		);

		ctx.estimate("distribution::gamma",
			[](real x) {return distribution::gamma(x, 3, 2);},
			[](real x) {return 4 * th::pow(x, 2) * th::exp(-2 * x);},
			exp_opt
		);

		ctx.estimate("distribution::beta",
			[](real x) {return distribution::beta(x, 2, 3);},
			[](real x) {return 12 * x * th::square(1 - x);},
			beta_opt
		);

		ctx.estimate("distribution::exponential",
			[](real x) {return distribution::exponential(x, 2);},
			[](real x) {return 2 * th::exp(-2 * x);},
			exp_opt
		);


		// Evaluate PDFs at known values

		ctx.equals("distribution::gaussian", distribution::gaussian(0, 0, 1), 1 / SQRT2 / SQRTPI);
		ctx.equals("distribution::gaussian (wrapper)", distribution::gaussian(0, vec<real>{0, 1}), 1 / SQRT2 / SQRTPI);

		ctx.equals("distribution::bernoulli", distribution::bernoulli(1, 0.25), 0.25);
		ctx.equals("distribution::bernoulli (wrapper)", distribution::bernoulli(0, vec<real>{0.25}), 0.75);

		ctx.equals("distribution::poisson", distribution::poisson(2, 3), 9 * th::exp(-3.0) / 2);
		ctx.equals("distribution::poisson (wrapper)", distribution::poisson(2, vec<real>{3}), 9 * th::exp(-3.0) / 2);

		ctx.equals("distribution::binomial", distribution::binomial(2, 4, 0.5), 0.375);
		ctx.equals("distribution::binomial (wrapper)", distribution::binomial(2, vec<real>{4, 0.5}), 0.375);

		ctx.equals("distribution::multinomial", distribution::multinomial({1, 2}, 3, 2, {0.25, 0.75}), 27.0 / 64.0);

		ctx.equals("distribution::chi_squared", distribution::chi_squared(1, 2), th::exp(-0.5) / 2);
		ctx.equals("distribution::chi_squared", distribution::chi_squared(1, 2, special::half_gamma(2)), th::exp(-0.5) / 2);
		ctx.equals("distribution::chi_squared (wrapper)", distribution::chi_squared(1, vec<real>{2}), th::exp(-0.5) / 2);

		ctx.equals("distribution::gamma (wrapper)", distribution::gamma(1, vec<real>{3, 2}), distribution::gamma(1, 3, 2));
		ctx.equals("distribution::beta (wrapper)", distribution::beta(0.5, vec<real>{2, 3}), distribution::beta(0.5, 2, 3));

		ctx.equals("distribution::student", distribution::student(0, 5), 0.3796066898, 1E-8);
		ctx.equals("distribution::student (wrapper)", distribution::student(0, vec<real>{5}), distribution::student(0, 5));

		ctx.equals("distribution::log_normal", distribution::log_normal(1, 0, 1), 1 / SQRT2 / SQRTPI);
		ctx.equals("distribution::log_normal (wrapper)", distribution::log_normal(1, vec<real>{0, 1}), distribution::log_normal(1, 0, 1));

		ctx.equals("distribution::rayleigh", distribution::rayleigh(1, 2), th::exp(-0.125) / 4);
		ctx.equals("distribution::rayleigh (wrapper)", distribution::rayleigh(1, vec<real>{2}), distribution::rayleigh(1, 2));

		ctx.equals("distribution::cauchy", distribution::cauchy(1, 0, 2), 1.0 / (PI * 2.5));
		ctx.equals("distribution::cauchy (wrapper)", distribution::cauchy(1, vec<real>{0, 2}), distribution::cauchy(1, 0, 2));

		ctx.equals("distribution::breit_wigner", distribution::breit_wigner(1, 0, 2), 1.0 / (2 * PI));
		ctx.equals("distribution::breit_wigner (wrapper)", distribution::breit_wigner(1, vec<real>{0, 2}), distribution::breit_wigner(1, 0, 2));

		ctx.equals("distribution::maxwell", distribution::maxwell(1, 2), (SQRT2 / SQRTPI) * th::exp(-0.125) / 8);
		ctx.equals("distribution::maxwell (wrapper)", distribution::maxwell(1, vec<real>{2}), distribution::maxwell(1, 2));

		ctx.equals("distribution::laplace", distribution::laplace(1, 0, 2), th::exp(-0.5) / 4);
		ctx.equals("distribution::laplace (wrapper)", distribution::laplace(1, vec<real>{0, 2}), distribution::laplace(1, 0, 2));

		ctx.equals("distribution::pareto", distribution::pareto(2, 1, 3), 3.0 / 16.0);
		ctx.equals("distribution::pareto (wrapper)", distribution::pareto(2, vec<real>{1, 3}), distribution::pareto(2, 1, 3));

		ctx.equals("distribution::erlang", distribution::erlang(2, 3, 1), 2 * th::exp(-2.0));
		ctx.equals("distribution::erlang (wrapper)", distribution::erlang(2, vec<real>{3, 1}), distribution::erlang(2, 3, 1));
	}


	// histogram.h
	{
		histogram hist (2, 0, 4);

		hist.insert(1);
		hist.insert(3);
		ctx.equals("histogram::index", hist.index(4), 1);

		ctx.equals(
			"histogram::range",
			(hist.range() - vec2({0, 4})).norm(),
			0.0
		);

		ctx.equals("histogram::number", hist.number(), 2);
		ctx.equals("histogram::bins[0]", hist.bins()[0], 1);
		ctx.equals("histogram::bins[1]", hist.bins()[1], 1);
		ctx.equals("histogram::range_lower", hist.range_lower(), 0.0);
		ctx.equals("histogram::range_upper", hist.range_upper(), 4.0);
		ctx.equals("histogram::max", hist.max(), 3.0);
		ctx.equals("histogram::min", hist.min(), 1.0);
		ctx.equals("histogram::mean", hist.mean(), 2.0);
		ctx.equals("histogram::tss", hist.tss(), 2.0);
		ctx.equals("histogram::operator()", hist(1), 1);
		ctx.equals("histogram::operator[]", hist[1], 1);
		ctx.equals("histogram::to_string", hist.to_string().empty(), false);

		ctx.equals("stats::mean(histogram)", stats::mean(hist), 2.0);
		ctx.equals("stats::tss(histogram)", stats::tss(hist), 2.0);
		ctx.equals("stats::variance(histogram)", stats::variance(hist), 2.0);
		ctx.equals("stats::stdev(histogram)", stats::stdev(hist), th::sqrt(2.0));

		ctx.equals("max(histogram)", theoretica::max(hist), 3.0);
		ctx.equals("min(histogram)", theoretica::min(hist), 1.0);

		hist.rebuild({2, 1}, vec2({0, 4}), 3, 2.0, 2.0, 1.0, 3.0);
		ctx.equals("histogram::rebuild", hist.number(), 3);

		// Test histogram construction from data
		histogram data_histogram (hist_data);
		ctx.equals("histogram(data)::number", data_histogram.number(), 12);
		ctx.equals("histogram(data)::mean", data_histogram.mean(), 431.0 / 12.0);
		ctx.equals("histogram(data)::tss", data_histogram.tss(), 5507.0 / 12.0);
		ctx.equals("histogram(data)::range_lower", data_histogram.range_lower(), 28.0);
		ctx.equals("histogram(data)::range_upper", data_histogram.range_upper(), 48.0);
		ctx.equals("histogram(data)::bins[0]", data_histogram[0], 6);
		ctx.equals("histogram(data)::bins[1]", data_histogram[1], 3);
		ctx.equals("histogram(data)::bins[2]", data_histogram[2], 3);


		// Test default range of uninitialized histogram (should contain NaN)
		histogram empty_histogram;
		ctx.equals("histogram().range()[0]", std::isnan(empty_histogram.range()[0]), true);
		ctx.equals("histogram().range()[1]", std::isnan(empty_histogram.range()[1]), true);

		// Test insertion behavior for values outside histogram range,
		// the histogram should silently reject out-of-range values
		const int count_before = hist.number();
		hist.insert(-1);
		hist.insert(5);
		const int count_after = hist.number();
		ctx.equals("histogram out-of-range", count_after, count_before);

		std::stringstream stream;
		stream << hist;
		ctx.equals("histogram::operator<<", stream.str().empty(), false);
		ctx.equals("histogram::to_string", hist.to_string().empty(), false);
	}


	// errorprop.h
	{
		auto product = [](vec<multidual<2>> values) {
			return values[0] * values[1];
		};

		const vec2 best = {2, 3};
		const vec2 errors = {0.1, 0.2};
		ctx.equals("stats::propagerr (independent)",
			stats::propagerr<2>(product, best, errors),
			0.5
		);

		ctx.equals(
			"stats::covar_mat",
			stats::covar_mat<mat<real, 2, 2>>(std::vector<vec<real>>{data, data2})(0, 1),
			covariance1
		);

		ctx.equals("stats::propagerr (covariance)",
			stats::propagerr<2>(
				product, best,
				mat<real, 2, 2>(vec<real>{0.01, 0, 0, 0.04})
			),
			0.5);

		ctx.equals(
			"stats::propagerr (datasets)",
			stats::propagerr<2>(product, std::vector<vec2>{{2, 2}, {3, 3}}),
			0.0, 1E-8
		);


		// Monte Carlo error propagation: sum of two independent variables
		// variance = 0.01^2 + 0.01^2 = 0.0002
		PRNG generator = PRNG::xoshiro(42);

		std::vector<pdf_sampler> samplers {
			pdf_sampler::gaussian(1, 0.01, generator), 
			pdf_sampler::gaussian(2, 0.01, generator) 
		};

		auto sum = [](vec<real> values) {return values[0] + values[1];};

		ctx.equals("stats::propagerr_mc", 
			stats::propagerr_mc(sum, samplers, 1E+06),
			th::sqrt(2E-04), 1E-02
		);
	}
}
