
///
/// @file attractor.cpp Compute the orbit of a chaotic attractor
///

#include <fstream>
#include "theoretica.h"
using namespace th;


// You can use the following gnuplot command to plot the results:
// splot "attractor.dat" using 2:3:4 with lines title "Attractor"


// Output data filename
const std::string filename = "attractor.dat";

// Initial conditions of the system
const vec3 initial_conditions = {0.1, 0, 0};

// Timestep
const real timestep = 0.001;

// Number of iterations
const int iterations = 100000;

// System parameters
const real a = 13;
const real b = 20;
const real c = 8 / 3.0;

// Unused constants
// const real d = 0;
// const real e = 0;
// const real g = 0;


// System of differential equations
vec3 f(real t, vec3 v) {

    // x, y, z
    const real x = v[0];
    const real y = v[1];
    const real z = v[2];

    // dx/dt, dy/dy, dz/dt
    return {
        a * y - a * x,
        x * b - x * z,
        x * y - c * z
    };
}


int main() {
    
    // Initialize the state of the system
    ode_state<3> state(initial_conditions);

    // Output file
    std::ofstream file(filename);

    for (int i = 0; i < iterations; ++i) {
        
        // Update the state of the system
        state = ode_rk4(f, state, timestep);

        // Print the state to file
        file << state << std::endl;
    }

    return 0;
}
