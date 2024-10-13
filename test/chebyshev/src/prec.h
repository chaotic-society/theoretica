
///
/// @file prec.h Precision testing module.
///

#ifndef CHEBYSHEV_PREC_H
#define CHEBYSHEV_PREC_H

#include <string>
#include <vector>
#include <map>
#include <fstream>
#include <iostream>
#include <ctime>

#include "./prec/prec_structures.h"
#include "./prec/fail.h"
#include "./prec/estimator.h"
#include "./core/output.h"
#include "./core/random.h"


namespace chebyshev {

	/// @namespace chebyshev::prec Precision testing module.
	///
	/// This module provides functions to estimate the precision
	/// and accuracy of mathematical approximations, over an entire
	/// domain using prec::estimate or at single points using prec:equals.
	/// For estimates over a domain, precision estimators are used.
	namespace prec {


		/// @class prec_settings Global settings of the precision testing module.
		struct prec_settings {
			
			/// Name of the module being tested
			std::string moduleName = "unknown";
			
			/// Print to standard output or not
			bool quiet = false;

			/// Output to file?
			bool outputToFile = true;

			/// Default number of iterations for integral quadrature
			unsigned int defaultIterations = CHEBYSHEV_PREC_ITER;

			/// Default fail function
			FailFunction defaultFailFunction = fail::fail_on_max_err();

			/// Default tolerance on max absolute error
			long double defaultTolerance = CHEBYSHEV_PREC_TOLERANCE;

			/// The files to write all precision testing results to
			std::vector<std::string> outputFiles {};

			/// Default columns to print for precision estimates.
			std::vector<std::string> estimateColumns = {
				"name", "meanErr", "rmsErr", "maxErr", "failed"
			};

			/// The files to write estimate results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> estimateOutputFiles {};

			/// Default columns to print for equations.
			std::vector<std::string> equationColumns = {
				"name", "difference", "tolerance", "failed"
			};

			/// The files to write equation results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> equationOutputFiles {};

			/// Target tests marked for execution,
			/// can be picked by passing test case names
			/// by command line. (all tests will be executed if empty)
			std::map<std::string, bool> pickedTests {};

		} settings;


		/// @class prec_results Test results of the precision testing module.
		struct prec_results {

			/// Total number of tests run
			unsigned int totalTests = 0;

			/// Number of failed tests
			unsigned int failedTests = 0;

			/// Results of error estimation
			std::map<std::string, std::vector<estimate_result>> estimateResults {};

			/// Results of equation evaluation
			std::map<std::string, std::vector<equation_result>> equationResults {};
			
		} results;


		/// Setup the precision testing environment.
		///
		/// @param moduleName Name of the module under test.
		/// @param argc The number of command line arguments
		/// @param argv A list of C-style strings containing
		/// the command line arguments.
		inline void setup(
			std::string moduleName,
			int argc = 0,
			const char** argv = nullptr) {


			// Initialize list of picked tests
			if(argc && argv)
				for (int i = 1; i < argc; ++i)
					settings.pickedTests[argv[i]] = true;

			std::cout << "Starting precision testing of the "
				<< moduleName << " module ..." << std::endl;

			settings.moduleName = moduleName;
			results.failedTests = 0;
			results.totalTests = 0;

			random::setup();
			output::setup();
		}


		/// Terminate the precision testing environment,
		/// printing the results to standard output and output files.
		///
		/// @param exit Whether to exit after terminating the module.
		inline void terminate(bool exit = true) {

			output::settings.quiet = settings.quiet;

			// Output to file is true but no specific files are specified, add default output file.
			if(	 settings.outputToFile &&
				!output::settings.outputFiles.size() &&
				!settings.estimateOutputFiles.size() &&
				!settings.equationOutputFiles.size() &&
				!settings.outputFiles.size()) {

				settings.outputFiles = { settings.moduleName + "_results" };
			}

			std::vector<std::string> outputFiles;

			// Print estimate results
			outputFiles  = settings.outputFiles;
			outputFiles.insert(outputFiles.end(), settings.estimateOutputFiles.begin(), settings.estimateOutputFiles.end());

			output::print_results(results.estimateResults, settings.estimateColumns, outputFiles);

			// Print equation results
			outputFiles  = settings.outputFiles;
			outputFiles.insert(outputFiles.end(), settings.equationOutputFiles.begin(), settings.equationOutputFiles.end());

			output::print_results(results.equationResults, settings.equationColumns, outputFiles);

			std::cout << "Finished testing " << settings.moduleName << '\n';
			std::cout << results.totalTests << " total tests, "
				<< results.failedTests << " failed (" << std::setprecision(3) <<
				(results.failedTests / (double) results.totalTests) * 100 << "%)"
				<< '\n';

			// Discard previous results
			results = prec_results();

			if(exit) {
				output::terminate();
				std::exit(results.failedTests);
			}
		}


		/// Estimate error integrals over a function
		/// with respect to an exact function,
		/// with the given options.
		///
		/// @param name The name of the test case
		/// @param funcApprox The approximation to test
		/// @param funcExpected The expected result
		/// @param opt The options for the estimation
		template<typename R, typename ...Args,
			typename Function1 = std::function<R(Args...)>,
			typename Function2 = Function1>
			
		inline void estimate(
			const std::string& name,
			Function1 funcApprox,
			Function2 funcExpected,
			estimate_options<R, Args...> opt) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(settings.pickedTests.size())
				if(settings.pickedTests.find(name) == settings.pickedTests.end())
					return;

			// Use the estimator to estimate error integrals.
			auto res = opt.estimator(funcApprox, funcExpected, opt);

			res.name = name;
			res.domain = opt.domain;
			res.tolerance = opt.tolerance;
			res.quiet = opt.quiet;
			res.iterations = opt.iterations;

			// Use the fail function to determine whether the test failed.
			res.failed = opt.fail(res);

			results.totalTests++;
			if(res.failed)
				results.failedTests++;

			results.estimateResults[name].push_back(res);
		}


		/// Estimate error integrals over a function
		/// with respect to an exact function.
		///
		/// @param name The name of the test case.
		/// @param funcApprox The approximation to test.
		/// @param funcExpected The expected result.
		/// @param intervals The (potentially multidimensional)
		/// domain of estimation.
		/// @param iterations The number of function evaluations.
		/// @param fail The fail function to determine whether
		/// the test failed.
		/// @param estimator The precision estimator to use.
		/// @param quiet Whether to output the result.
		template<typename R, typename ...Args,
			typename Function1 = std::function<R(Args...)>,
			typename Function2 = Function1>

		inline void estimate(
			const std::string& name,
			Function1 funcApprox,
			Function2 funcExpected,
			std::vector<interval> domain,
			long double tolerance, unsigned int iterations,
			FailFunction fail,
			Estimator<R, Args...> estimator,
			bool quiet = false) {

			estimate_options<R, Args...> opt {};
			opt.domain = domain;
			opt.tolerance = tolerance;
			opt.iterations = iterations;
			opt.fail = fail;
			opt.estimator = estimator;
			opt.quiet = quiet;

			estimate(name, funcApprox, funcExpected, opt);
		}


		/// Estimate error integrals over a real function
		/// of real variable, with respect to an exact function.
		///
		/// @param name The name of the test case.
		/// @param funcApprox The approximation to test.
		/// @param funcExpected The expected result.
		/// @param intervals The (potentially multidimensional)
		/// domain of estimation.
		/// @param iterations The number of function evaluations.
		/// @param fail The fail function to determine whether
		/// the test failed (defaults to fail_on_max_err).
		/// @param estimator The precision estimator to use
		/// (defaults to the trapezoid<double> estimator).
		/// @param quiet Whether to output the result.
		inline void estimate(
			const std::string& name,
			EndoFunction<double> funcApprox,
			EndoFunction<double> funcExpected,
			interval domain,
			long double tolerance = settings.defaultTolerance,
			unsigned int iterations = settings.defaultIterations,
			FailFunction fail = fail::fail_on_max_err(),
			Estimator<double, double> estimator = estimator::quadrature1D<double>(),
			bool quiet = false) {

			estimate_options<double, double> opt {};
			opt.domain = { domain };
			opt.tolerance = tolerance;
			opt.iterations = iterations;
			opt.fail = fail;
			opt.estimator = estimator;
			opt.quiet = quiet;

			estimate(name, funcApprox, funcExpected, opt);
		}

		/// @namespace chebyshev::prec::property Property testing of functions
		///
		/// When estimating error integrals, it is usually necessary to have
		/// a function to compare the result to, considered exact. Using
		/// property testing, it is possible to test a specific property
		/// of a function (such as involution or homogeneity) doing away
		/// with the additional "exact" function.
		namespace property {

			/// Precision testing of an endofunction which is
			/// equivalent to the identity.
			///
			/// @param name The name of the test case.
			/// @param id The identity function to test.
			/// @param opt The options for estimation.
			template<typename Type, typename Identity = EndoFunction<Type>>
			inline void identity(
				const std::string& name,
				Identity id,
				const estimate_options<Type, Type>& opt) {

				// Apply the identity function
				EndoFunction<Type> funcApprox = [&](Type x) -> Type {
					return id(x);
				};

				// And compare it to the identity
				EndoFunction<Type> funcExpected = [](Type x) -> Type {
					return x;
				};

				estimate(name, funcApprox, funcExpected, opt);
			}

			/// Precision testing of an endofunction which is
			/// an involution. The function is applied two times
			/// to input values and it is checked against the identity.
			///
			/// @param name The name of the test case.
			/// @param involution The involution to test.
			/// @param opt The options for estimation.
			template<typename Type, typename Involution = EndoFunction<Type>>
			inline void involution(
				const std::string& name,
				Involution invol,
				const estimate_options<Type, Type>& opt) {

				// Apply the involution two times
				EndoFunction<Type> funcApprox = [&](Type x) -> Type {
					return invol(invol(x));
				};

				// And compare it to the identity
				EndoFunction<Type> funcExpected = [](Type x) -> Type {
					return x;
				};

				estimate(name, funcApprox, funcExpected, opt);
			}


			/// Precision testing of an endofunction which is
			/// idempotent. The function is applied two times
			/// to input values and it is checked against itself.
			///
			/// @param name The name of the test case.
			/// @param idem The idempotent function to test.
			/// @param opt The options for estimation.
			template<typename Type, typename Involution = EndoFunction<Type>>
			inline void idempotence(
				const std::string& name,
				Involution idem,
				const estimate_options<Type, Type>& opt) {

				// Apply the idem two times
				EndoFunction<Type> funcApprox = [&](Type x) -> Type {
					return idem(idem(x));
				};

				// And compare it to the identity
				EndoFunction<Type> funcExpected = [&](Type x) -> Type {
					return idem(x);
				};

				estimate(name, funcApprox, funcExpected, opt);
			}


			/// Precision testing of an function which is
			/// homogeneous over the domain. The function is applied
			/// to input values and it is checked against zero.
			/// The zero value is constructed as OutputType(0.0), but may
			/// be specified as an additional argument.
			///
			/// @param name The name of the test case.
			/// @param hom The homogeneous function to test.
			/// @param opt The options for estimation.
			/// @param zero_element The zero element of type OutputType
			/// (defaults to OutputType(0.0)).
			template<typename InputType, typename OutputType = InputType,
			typename Homogeneous = std::function<OutputType(InputType)>>
			inline void homogeneous(
				const std::string& name,
				Homogeneous hom,
				const estimate_options<InputType, OutputType>& opt,
				OutputType zero_element = OutputType(0.0)) {

				// Apply the homogeneous function
				std::function<OutputType(InputType)> funcApprox =
					[&](InputType x) -> OutputType {
						return hom(x);
					};

				// And compare it to the zero element
				std::function<OutputType(InputType)> funcExpected =
					[&](InputType x) -> OutputType {
						return zero_element;
					};

				estimate(name, funcApprox, funcExpected, opt);
			}
		}


		/// Test an equivalence up to a tolerance,
		/// with the given options (e.g. for residual testing).
		///
		/// @param name The name of the test case
		/// @param evaluate The evaluated value
		/// @param expected The expected value
		/// @param opt The options for the evaluation
		template<typename T = double>
		inline void equals(
			const std::string& name,
			const T& evaluated, const T& expected,
			equation_options<T> opt = equation_options<T>()) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(settings.pickedTests.size())
				if(settings.pickedTests.find(name) == settings.pickedTests.end())
					return;

			equation_result res {};

			long double diff = opt.distance(evaluated, expected);

			// Mark the test as failed if the
			// distance between the two values
			// is bigger than the tolerance.
			res.failed = (diff > opt.tolerance);

			res.name = name;
			res.difference = diff;
			res.tolerance = opt.tolerance;
			res.quiet = opt.quiet;

			results.totalTests++;
			if(res.failed)
				results.failedTests++;

			// Register the result of the equation by name
			results.equationResults[name].push_back(res);
		}


		/// Test an equivalence up to a tolerance,
		/// with the given options (e.g. for residual testing).
		///
		/// @param name The name of the test case
		/// @param evaluate The evaluated value
		/// @param expected The expected value
		/// @param distance The distance function to use
		/// @param tolerance The tolerance for the evaluation
		/// @param quiet Whether to output the result
		template<typename T = double>
		inline void equals(
			const std::string& name,
			const T& evaluated, const T& expected,
			long double tolerance,
			DistanceFunction<T> distance,
			bool quiet = false) {

			equation_options<T> opt {};
			opt.tolerance = tolerance;
			opt.distance = distance;
			opt.quiet = quiet;

			equals(name, evaluated, expected, opt);
		}


		/// Test an equivalence up to a tolerance,
		/// with the given options (e.g. for residual testing).
		///
		/// @param name The name of the test case
		/// @param evaluate The evaluated value
		/// @param expected The expected value
		/// @param tolerance The tolerance for the evaluation
		/// @param quiet Whether to output the result
		inline void equals(
			const std::string& name,
			long double evaluated, long double expected,
			long double tolerance = settings.defaultTolerance,
			bool quiet = false) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(settings.pickedTests.size())
				if(settings.pickedTests.find(name) == settings.pickedTests.end())
					return;

			equation_result res {};

			long double diff = distance::abs_distance(evaluated, expected);

			// Mark the test as failed if the
			// distance between the two values
			// is bigger than the tolerance.
			res.failed = (diff > tolerance);

			res.name = name;
			res.difference = diff;
			res.tolerance = tolerance;
			res.quiet = quiet;

			res.evaluated = evaluated;
			res.expected = expected;

			results.totalTests++;
			if(res.failed)
				results.failedTests++;

			// Register the result of the equation by name
			results.equationResults[name].push_back(res);
		}


		/// Evaluate multiple pairs of values for equivalence
		/// up to the given tolerance (e.g. for residual testing).
		///
		/// @param name The name of the function or test case
		/// @param values A list of values to equate
		/// @param tolerance The tolerance for the evaluation
		/// @param quiet Whether to output the result
		template<typename T>
		inline void equals(
			const std::string& name,
			std::vector<std::array<T, 2>> values,
			long double tolerance = settings.defaultTolerance,
			bool quiet = false) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(settings.pickedTests.size())
				if(settings.pickedTests.find(name) == settings.pickedTests.end())
					return;

			for (const auto& v : values)
				equals(name, v[0], v[1], tolerance, quiet);
		}
	}
}

#endif
