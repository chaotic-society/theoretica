
///
/// @file attractor.cpp Compute the orbit of a chaotic attractor
///

#include <fstream>
#include <iostream>

#include "theoretica.h"
using namespace th;


// You can use the following gnuplot command to plot the results:
// splot "attractor.dat" using 2:3:4 with lines title "Attractor"


// Output data filename
const std::string filename = "examples/attractor.dat";

// Initial conditions of the system
const vec3 x0 = {0.1, 0, 0};

// Starting time
const real t0 = 0.0;

// Final time
const real tf = 50.0;

// Timestep
const real timestep = 0.001;

// System parameters
const real a = 13;
const real b = 20;
const real c = 8./3.;


// System of differential equations
vec3 f(real t, vec3 v) {

    const real x = v[0];
    const real y = v[1];
    const real z = v[2];

    // Vector of time derivatives
    return {
        a * y - a * x,
        x * b - x * z,
        x * y - c * z
    };
}


int main() {

    // Solve the system of differential equations using Runge-Kutta's method
    auto solution = ode::solve_rk4(
        f, x0, t0, tf, timestep
    );

    std::ofstream file (filename);
    
    if (!file.is_open()) {
        std::cout << "Unable to open file: " << filename << std::endl;
        return 1;
    }

    // Write the solution to file
    file << solution;
}
