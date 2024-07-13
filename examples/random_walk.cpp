
///
/// @file random_walk.cpp A 3D random walk generator
/// This example may be compiled using 'make random_walk'
///

#include "theoretica.h"
using namespace th;

#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>


// You can plot the data file using gnuplot:
// splot "examples/random_walk.dat" with lines title "3D Random Walk"


int main() {

    // Number of steps
    const unsigned int N = 10000;

    // Output file
    std::ofstream file("examples/random_walk.dat");

    // Random number generator
    PRNG g = PRNG::xoshiro(time(nullptr));

    // Trajectory
    std::vector<vec3> pos(N);

    // Start from the origin
    pos[0] = {0, 0, 0};
    file << pos[0].to_string(" ", false) << std::endl;

    for (unsigned int i = 1; i < N; ++i) {

        // Generate a random vector with fixed length
        // You can change the generation of r to test
        // different distributions (e.g. rand_gaussian(0, 1, g))
        const real r = 1;
        const real theta = rand_uniform(0, PI, g);
        const real phi = rand_uniform(0, TAU, g);

        // Update the trajectory
        pos[i] = pos[i - 1] + vec3({
            r * th::sin(theta) * th::cos(phi),
            r * th::sin(theta) * th::sin(phi),
            r * th::cos(theta)
        });

        file << pos[i].to_string(" ", false) << std::endl;
    }


    // Compute the mean direction of the steps
    vec3 mean_dir = {0, 0, 0};
    for (unsigned int i = 1; i < N; ++i)
        mean_dir += (pos[i] - pos[i - 1]) / N;

    std::cout << "Mean Direction: " << mean_dir << std::endl;


    // Compute the RMS displacement
    real sqr_disp = 0;
    for (unsigned int i = 0; i < N; ++i)
        sqr_disp += pos[i] * pos[i] / N;

    std::cout << "RMS Displacement: " << th::sqrt(sqr_disp) << std::endl;

    return 0;
}
