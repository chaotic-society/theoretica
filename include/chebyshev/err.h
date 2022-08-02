
///
/// @file error_handling.h Error handling checks
///

#pragma once

#include "core/common.h"
#include <vector>
#include <cstdlib>


#define AUTOASSERT(exp) chebyshev::err::assert(exp, #exp);


namespace chebyshev {

	/// @namespace chebyshev::err Error testing
	namespace err {


		/// @class err_state Global state of error testing
		struct err_state {
				
			/// Name of the module being tested
			std::string moduleName = "unknown";

			/// Print to standard output?
			bool quiet = false;

			// Total number of checks
			unsigned int totalChecks = 0;

			/// Number of failed checks
			unsigned int failedChecks = 0;

		} state;


		/// Setup error checking
		void setup(const std::string& module) {

			state.moduleName = module;
			srand(time(nullptr));

			std::cout << "Starting error checking on " << state.moduleName << " ...\n" << std::endl;
		}


		/// Terminate error checking
		void terminate(bool exit = true) {

			std::cout << "\nEnding error checking on " << state.moduleName << " ..." << std::endl;
			std::cout << state.totalChecks << " total checks, " << state.failedChecks << " failed ("
				<< (state.failedChecks / Real(state.totalChecks) * 100) << "%)" << std::endl;

			if(exit)
				std::exit(state.failedChecks);
		}


		/// Assert that an expression is true
		void assert(bool exp, std::string desc = "") {

			state.totalChecks++;

			if(!exp) {
				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks << ")";

					if(desc != "")
						std::cout << ":\n\t\tDescription: " << desc << std::endl;
					else
						std::cout << " (no description provided)" << std::endl;

					std::cout << std::endl;
				}
				
			} else if(!state.quiet) {
				std::cout << "\tSuccessful assert (n. " << state.totalChecks << ")";
				if(desc != "")
					std::cout << ":\n\t\tDescription: " << desc << std::endl;
				else
					std::cout << " (no description provided)" << std::endl;

				std::cout << std::endl;
			}
		}


		/// Check errno value after function call
		void check_errno(RealFunction f, Real x, int exp_errno) {

			state.totalChecks++;

			try {
				volatile Real r = f(x);
				r += 1;
			} catch(...) {}

			if(errno != exp_errno) {

				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks << ")" << std::endl;
					std::cout << "\t\tExpected ERRNO: " << exp_errno << std::endl;
					std::cout << "\t\tEvaluated ERRNO: " << errno << std::endl;
					std::cout << "\t\tInput: " << x << std::endl;
					std::cout << std::endl;
				}
				
			} else {

				if(!state.quiet) {
					std::cout << "\tSuccessful assert (n. " << state.totalChecks
						<< "): ERRNO was set correctly" << std::endl;
				}

			}

		}


		/// Check errno value after function call
		void check_errno(RealFunction f, RealInputGenerator g, int exp_errno) {

			state.totalChecks++;

			Real x = g(rand());

			try {
				volatile Real r = f(x);
				r += 1;
			} catch(...) {}

			if(errno != exp_errno) {

				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks << ")" << std::endl;
					std::cout << "\t\tExpected ERRNO: " << exp_errno << std::endl;
					std::cout << "\t\tEvaluated ERRNO: " << errno << std::endl;
					std::cout << std::endl;
				}
				
			} else {

				if(!state.quiet) {
					std::cout << "\tSuccessful assert (n. " << state.totalChecks
						<< "): ERRNO was set correctly" << std::endl;
				}

			}

		}


		/// Check errno value after function call
		void check_errno(RealFunction f, RealInputGenerator g, const std::vector<int>& exp_flags) {


			state.totalChecks++;

			Real x = g(rand());
			bool all_set = true;

			try {
				volatile Real r = f(x);
				r += 1;
			} catch(...) {}

			for (unsigned int i = 0; i < exp_flags.size(); ++i) {
				if(!(errno & exp_flags[i]))
					all_set = false;
			}

			if(!all_set) {

				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks
						<< "): ERRNO was NOT set correctly" << std::endl;
					std::cout << std::endl;
				}
				
			} else {

				if(!state.quiet) {
					std::cout << "\tSuccessful assert (n. " << state.totalChecks
						<< "): ERRNO was set correctly" << std::endl;
				}

			}

		}


		/// Check that an exception is thrown during a function call
		void check_exception(RealFunction f, Real x) {

			state.totalChecks++;

			bool thrown = false;

			try {
				volatile Real r = f(x);
				r += 1;
			} catch(...) {
				thrown = true;
			}

			if(!thrown) {

				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks << ")" << std::endl;
					std::cout << "\tNo exception was thrown" << std::endl;
					std::cout << std::endl;
				}
				
			} else if(!state.quiet) {
				std::cout << "\tSuccessful assert (n. " << state.totalChecks
					<< "): Exception was thrown correctly" << std::endl;
			}

		}


		/// Check that an exception is thrown during a function call
		void check_exception(RealFunction f, RealInputGenerator g) {

			state.totalChecks++;

			bool thrown = false;
			Real x = g(rand());

			try {
				volatile Real r = f(x);
				r += 1;
			} catch(...) {
				thrown = true;
			}

			if(!thrown) {

				state.failedChecks++;

				if(!state.quiet) {
					std::cout << "\tFailed assert (n. " << state.totalChecks << "):" << std::endl;
					std::cout << "\tNo exception was thrown" << std::endl;
					std::cout << std::endl;
				}
				
			} else if(!state.quiet) {
				std::cout << "\tSuccessful assert (n. " << state.totalChecks
					<< "): Exception was thrown correctly" << std::endl;
			}

		}

	}

}
