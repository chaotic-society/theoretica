
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
	namespace prec {


		/// @class prec_state
		/// Global state of the precision testing module.
		struct prec_state {
			
			/// Name of the module being tested
			std::string moduleName = "unknown";
			
			/// Print to standard output or not
			bool quiet = false;

			/// Output to file?
			bool outputToFile = true;

			/// Total number of tests run
			unsigned int totalTests = 0;

			/// Number of failed tests
			unsigned int failedTests = 0;

			/// Default number of iterations for integral quadrature
			unsigned int defaultIterations = CHEBYSHEV_PREC_ITER;

			/// Default fail function
			FailFunction defaultFailFunction = fail::fail_on_max_err();

			/// Default tolerance on max absolute error
			long double defaultTolerance = CHEBYSHEV_PREC_TOLERANCE;

			/// Results of precision testing
			std::map<std::string, std::vector<estimate_result>> estimateResults {};

			/// Default columns to print for precision estimates.
			std::vector<std::string> estimateColumns = {
				"funcName", "meanErr", "rmsErr", "maxErr", "failed"
			};

			/// The files to write estimate results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> estimateFiles {};
			
			/// Results of equations
			std::map<std::string, std::vector<equation_result>> equationResults {};

			/// Default columns to print for equations.
			std::vector<std::string> equationColumns = {
				"funcName", "difference", "tolerance", "failed"
			};

			/// The files to write equation results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> equationFiles {};

			/// Target tests marked for execution,
			/// can be picked by passing test case names
			/// by command line. (all tests will be executed if empty)
			std::map<std::string, bool> pickedTests {};

		} state;


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
					state.pickedTests[argv[i]] = true;

			std::cout << "Starting precision testing of the "
				<< moduleName << " module ..." << std::endl;

			state.moduleName = moduleName;
			state.failedTests = 0;
			state.totalTests = 0;

			random::setup();
			output::setup();
		}


		/// Terminate the precision testing environment.
		///
		/// @param exit Whether to exit after terminating the module.
		inline void terminate(bool exit = true) {

			output::state.quiet = state.quiet;

			// Output to file is true but no specific files are specified,
			// add default output file.
			if(state.outputToFile && !state.estimateFiles.size() && !state.equationFiles.size()) {
				std::string filename;
				filename = output::state.outputFolder + state.moduleName + "_results";
				output::state.outputFiles[filename] = std::ofstream(filename);
			}

			output::print_results(state.estimateResults, state.estimateColumns, state.estimateFiles);
			output::print_results(state.equationResults, state.equationColumns, state.equationFiles);

			std::cout << "Finished testing " << state.moduleName << '\n';
			std::cout << state.totalTests << " total tests, "
				<< state.failedTests << " failed (" <<
				(state.failedTests / (double) state.totalTests) * 100 << "%)"
				<< '\n';

			state = prec_state();

			if(exit) {
				output::terminate();
				std::exit(state.failedTests);
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
			if(state.pickedTests.size())
				if(state.pickedTests.find(name) == state.pickedTests.end())
					return;

			// Use the estimator to estimate error integrals.
			auto res = opt.estimator(funcApprox, funcExpected, opt);

			res.funcName = name;
			res.domain = opt.domain;
			res.tolerance = opt.tolerance;
			res.quiet = opt.quiet;

			// Use the fail function to determine whether the test failed.
			res.failed = opt.fail(res);

			state.totalTests++;
			if(res.failed)
				state.failedTests++;

			state.estimateResults[name].push_back(res);
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
			long double tolerance = CHEBYSHEV_PREC_TOLERANCE,
			unsigned int iterations = state.defaultIterations,
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

		namespace property {

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
			equation_options<T> opt) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(state.pickedTests.size())
				if(state.pickedTests.find(name) == state.pickedTests.end())
					return;

			equation_result res {};

			long double diff = opt.distance(evaluated, expected);

			// Mark the test as failed if the
			// distance between the two values
			// is bigger than the tolerance.
			res.failed = (diff > opt.tolerance);

			res.funcName = name;
			res.difference = diff;
			res.tolerance = opt.tolerance;
			res.quiet = opt.quiet;

			state.totalTests++;
			if(res.failed)
				state.failedTests++;

			// Register the result of the equation by name
			state.equationResults[name].push_back(res);
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
			long double tolerance = state.defaultTolerance,
			bool quiet = false) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(state.pickedTests.size())
				if(state.pickedTests.find(name) == state.pickedTests.end())
					return;

			equation_result res {};

			long double diff = distance::abs_distance(evaluated, expected);

			// Mark the test as failed if the
			// distance between the two values
			// is bigger than the tolerance.
			res.failed = (diff > tolerance);

			res.funcName = name;
			res.difference = diff;
			res.tolerance = tolerance;
			res.quiet = quiet;

			res.evaluated = evaluated;
			res.expected = expected;

			state.totalTests++;
			if(res.failed)
				state.failedTests++;

			// Register the result of the equation by name
			state.equationResults[name].push_back(res);
		}


		/// Evaluate multiple pairs of values for equivalence
		/// up to the given tolerance (e.g. for residual testing).
		inline void equals(
			const std::string& name,
			std::vector<std::array<long double, 2>> values,
			long double tolerance = state.defaultTolerance,
			bool quiet = false) {

			// Skip the test case if any tests have been picked
			// and this one was not picked.
			if(state.pickedTests.size())
				if(state.pickedTests.find(name) == state.pickedTests.end())
					return;

			for (const auto& v : values)
				equals(name, v[0], v[1], tolerance, quiet);
		}


	}
}

#endif
