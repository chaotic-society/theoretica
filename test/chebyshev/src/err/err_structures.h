
///
/// @file err_structures.h Structures for the error checking module.
///

#ifndef CHEBYSHEV_ERR_STRUCTURES_H
#define CHEBYSHEV_ERR_STRUCTURES_H

#include <string>
#include <vector>
#include <map>


namespace chebyshev {

	namespace err {


		/// @class assert_result
		/// Result of assertion checking of a function
		struct assert_result {

			/// Identifying name of the function or test case.
			std::string name = "unknown";
			
			/// Evaluated boolean value
			bool evaluated = false;

			/// Description of the assertion
			std::string description = "";

			/// Whether the test failed
			bool failed = true;

			/// Whether to print to standard output or not.
			bool quiet = false;
		};


		/// @class errno_result
		/// Result of errno checking of a function
		struct errno_result {

			/// Identifying name of the function or test case.
			std::string name = "unknown";
			
			/// Evaluated errno value
			int evaluated;

			/// Expected errno flags
			std::vector<int> expectedFlags;

			/// Whether the test failed.
			bool failed = true;

			/// Whether to print to standard output or not.
			bool quiet = false;
		};


		/// @class exception_result
		/// Result of exception checking of a function
		struct exception_result {

			/// Identifying name of the function or test case.
			std::string name = "unknown";
			
			/// Whether the function has thrown
			bool thrown = false;

			/// Whether the type of the thrown exception
			/// was correct.
			bool correctType = true;

			/// Whether the test failed.
			bool failed = true;

			/// Whether to print to standard output or not.
			bool quiet = false;
		};

	}
}

#endif
