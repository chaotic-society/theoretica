
///
/// @file montecarlo_integral.cpp A comparison between Monte Carlo techniques.
/// This example may be compiled using 'make montecarlo_integral'
///

#include <iostream>
#include "theoretica.h"
using namespace th;

#include <ctime>

// A simple function to integrate
real f(real x) {
    return th::sin(x);
}

vec2 F(vec2 x) {
    return {x[1] * th::sin(x[0]), x[0] * th::cos(x[1])};
}

int main() {


    // Initialize the pseudorandom number generator to be used
    random::XoshiroPrng g (time(nullptr));

    // Starting number of iterations
    unsigned int N = 10;

    // For pretty printing the results
    data_table table = data_table(
        0, {"N", "MC Error", "QMC Error"}
    );

    // Repeat the integration for different N (exponentially growing)
    for (; N <= 10'000'000; N *= 10) {

        // Crude Monte-Carlo
        auto crude = random::integral_mc(f, 0, PI/2, g, N);

        // Crude Quasi-Monte Carlo
        auto q_crude = random::integral_qmc(f, 0, PI/2, g, N);

        // The result, by default, is a stoch_result<real> object,
        // containing the estimate <value> and its <stdev>.
        // You can cast it to real, discarding additional information,
        // but this is not generally recommended for robust code.

        table["N"].append(N);
        table["MC Error"].append(crude.stdev);
        table["QMC Error"].append(q_crude.stdev);
    }

    io::println(table);

    // Integrate a multivariate function
    vec<vec2> domain = {{0, PI/2}, {0, PI/2}};
    io::println(random::integral_mc(F, domain, g));
}
