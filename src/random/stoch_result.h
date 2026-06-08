
///
/// @file stoch_result.h A structure holding information for stochastic methods.
///

#ifndef THEORETICA_STOCH_RESULT_H
#define THEORETICA_STOCH_RESULT_H

#ifndef THEORETICA_NO_PRINT
#include <sstream>
#include <ostream>
#endif

#include "../core/constants.h"
#include "../core/real_analysis.h"
#include "../core/core_traits.h"


namespace theoretica {

    /// @class stoch_result
    /// A structure holding information about the
    /// results of a stochastic method.
    ///
    /// @tparam Type The type of the computed result.
    template<typename Type>
    struct stoch_result {

        /// Estimated value
        Type value;

        /// Standard deviation of the estimate
        Type stdev;

        /// Sample size (0 if unused)
        size_t size {0};


        /// Automatically convert to Type
        /// to discard additional information.
        operator Type () const {
            return value;
        }


#ifndef THEORETICA_NO_PRINT

		/// Convert the stochastic result to string representation (scalar case)
        template<typename T = Type, disable_vector<T> = true>
        inline std::string to_string() const {

			std::stringstream res;
            res << value << " ± " << stdev << " (n = " << size << ")";
			return res.str();
		}

        /// Convert the stochastic result to string representation (vector case)
        template<typename T = Type, enable_vector<T> = true>
        inline std::string to_string() const {

			std::stringstream res;
            res << "(";

            for (unsigned int i = 0; i < value.size(); ++i)
                res << value[i] << " ± " << stdev[i] << (i < value.size() - 1 ? ", " : "");

            res << ")  (n = " << size << ")";
			return res.str();
		}

		/// Convert the stochastic result to string representation.
		inline operator std::string() {
			return to_string();
		}

		/// Stream the stochastic result in string representation
		/// to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const stoch_result& obj) {
			return out << obj.to_string();
		}

#endif

    };

}

#endif
