
///
/// @file err.h Error checking module
///

#ifndef CHEBYSHEV_ERR_H
#define CHEBYSHEV_ERR_H

#include <vector>
#include <cstdlib>
#include <iostream>

#include "./core/common.h"
#include "./core/random.h"
#include "./err/err_structures.h"


namespace chebyshev {


	// To err is human; to forgive, divine.


	/// @namespace chebyshev::err Error checking module
	namespace err {


		/// @class err_state
		/// Global state of the error testing module.
		struct err_state {
				
			/// Name of the module being tested
			std::string moduleName = "unknown";

			// Total number of checks
			unsigned int totalChecks = 0;

			/// Number of failed checks
			unsigned int failedChecks = 0;

			/// Whether to print to an output file.
			bool outputToFile = true;

			/// Results of checking assertions
			std::map<std::string, std::vector<assert_result>> assertResults {};

			/// The files to write assertion results results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> assertFiles {};

			/// Default columns to print for assertions.
			std::vector<std::string> assertColumns = {
				"funcName", "evaluated", "failed", "description"
			};

			/// Results of checking errno
			std::map<std::string, std::vector<errno_result>> errnoResults {};

			/// The files to write errno results results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> errnoFiles {};

			/// Default columns to print for errno checks.
			std::vector<std::string> errnoColumns = {
				"funcName", "evaluated", "expectedFlags", "failed"
			};

			/// Results of exception testing
			std::map<std::string, std::vector<exception_result>> exceptionResults {};

			/// The files to write exception results results to
			/// (if empty, all results are output to a generic file).
			std::vector<std::string> exceptionFiles {};

			/// Default columns to print for exception checks.
			std::vector<std::string> exceptionColumns = {
				"funcName", "thrown", "correctType", "failed"
			};

			/// Target checks marked for execution,
			/// can be picked by passing test case names
			/// by command line. (all tests will be executed if empty)
			std::map<std::string, bool> pickedChecks {};

			/// Whether to print to standard output
			bool quiet = false;

		} state;


		/// Setup error checking module.
		void setup(
			const std::string& moduleName, int argc = 0, const char** argv = nullptr) {

			if(argc && argv)
				for (int i = 1; i < argc; ++i)
					state.pickedChecks[argv[i]] = true;

			std::cout << "Starting error checking on "
				<< moduleName << " ..." << std::endl;

			state.moduleName = moduleName;
			state.failedChecks = 0;
			state.totalChecks = 0;

			random::setup();
			output::setup();
		}


		/// Terminate the error testing environment.
		///
		/// @param exit Whether to exit after terminating the module.
		void terminate(bool exit = true) {

			output::state.quiet = state.quiet;

			// Output to file is true but no specific files are specified,
			// add default output file.
			if(state.outputToFile && !state.assertFiles.size()
				&& !state.errnoFiles.size() && !state.exceptionFiles.size()) {
				std::string filename;
				filename = output::state.outputFolder + state.moduleName + "_results";
				output::state.outputFiles[filename] = std::ofstream(filename);
			}

			output::print_results(state.assertResults, state.assertColumns, state.assertFiles);
			output::print_results(state.errnoResults, state.errnoColumns, state.errnoFiles);
			output::print_results(state.exceptionResults, state.exceptionColumns, state.exceptionFiles);

			std::cout << "Finished error checking " << state.moduleName << " ...\n";
			std::cout << state.totalChecks
				<< " total checks, "
				<< state.failedChecks << " failed ("  << std::setprecision(3)
				<< (state.failedChecks / (double) state.totalChecks * 100.0)
				<< "%)" << std::endl;

			state = err_state();

			if(exit) {
				output::terminate();
				std::exit(state.failedChecks);
			}
		}


		/// Assert that an expression is true
		///
		/// @param name Name of the check (function name or test case name).
		/// @param exp Expression to test for truth.
		/// @param description Description of the assertion.
		void assert(const std::string& name, bool exp, std::string description = "") {

			assert_result res {};

			res.funcName = name;
			res.evaluated = exp;
			res.failed = !exp;
			res.description = description;

			state.totalChecks++;

			if(!exp)
				state.failedChecks++;

			state.assertResults[name].push_back(res);
		}


		/// Check errno value after function call
		template<typename Function, typename InputType>
		void check_errno(
			const std::string& name, Function f,
			InputType x, int expected_errno) {

			errno_result res {};
			errno = 0;

			try {
				volatile auto r = f(x);
				r = *(&r);
			} catch(...) {}

			res.funcName = name;
			res.evaluated = errno;
			res.expectedFlags = { expected_errno };
			res.failed = (errno != expected_errno);


			state.totalChecks++;

			if(res.failed)
				state.failedChecks++;

			state.errnoResults[name].push_back(res);
		}


		/// Check errno value after function call
		template<typename Function, typename InputType>
		void check_errno(
			const std::string& name, Function f,
			std::function<InputType(unsigned int)> generator,
			int expected_errno) {

			check_errno(name, f, generator(rand()), expected_errno);
		}


		/// Check errno value after function call
		template<typename Function, typename InputType>
		void check_errno(
			const std::string& name, Function f,
			InputType x, std::vector<int>& expected_flags) {


			errno_result res {};
			errno = 0;

			try {
				volatile auto r = f(x);
				r = *(&r);
			} catch(...) {}

			res.funcName = name;
			res.evaluated = errno;
			res.expectedFlags = expected_flags;
			
			res.failed = false;
			for (int flag : expected_flags)
				if(!(errno & flag))
					res.failed = true;

			state.totalChecks++;

			if(res.failed)
				state.failedChecks++;

			state.errnoResults[name].push_back(res);
		}


		/// Check errno value after function call
		template<typename Function, typename InputType>
		void check_errno(
			const std::string& name, Function f,
			std::function<InputType(unsigned int)> generator,
			std::vector<int>& expected_flags) {

			check_errno(name, f, generator(rand()), expected_flags);
		}


		/// Check that an exception is thrown during a function call
		template<typename Function, typename InputType>
		void check_exception(const std::string& name, Function f, InputType x) {

			exception_result res {};
			bool thrown = false;

			try {
				volatile auto r = f(x);
				r = *(&r);
			} catch(...) {
				thrown = true;
			}

			res.funcName = name;
			res.thrown = thrown;
			res.failed = !thrown;
			res.correctType = true;

			state.totalChecks++;
			if(!thrown)
				state.failedChecks++;

			state.exceptionResults[name].push_back(res);
		}


		/// Check that an exception is thrown during a function call
		template<typename Function, typename InputType>
		void check_exception(
			const std::string& name, Function f,
			std::function<InputType(unsigned int)> generator) {

			check_exception(name, f, generator(rand()));
		}


		template<typename ExceptionType, typename Function, typename InputType>
		void check_exception(const std::string& name, Function f, InputType x) {

			exception_result res {};
			bool thrown = false;
			bool correctType = false;

			try {
				volatile auto r = f(x);
				r = *(&r);
			} catch(ExceptionType& exc) {

				correctType = true;
				thrown = true;

			} catch(...) {
				thrown = true;
			}

			res.funcName = name;
			res.thrown = thrown;
			res.failed = !(thrown && correctType);
			res.correctType = correctType;

			state.totalChecks++;
			if(!thrown)
				state.failedChecks++;

			state.exceptionResults[name].push_back(res);
		}


		template<typename ExceptionType, typename Function, typename InputType>
		void check_exception(
			const std::string& name, Function f,
			std::function<InputType(unsigned int)> generator) {

			check_exception(name, f, generator(rand()));
		}

	}

}

#endif
