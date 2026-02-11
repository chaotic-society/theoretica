
///
/// @file histogram.cpp Read data points from the given file and construct an histogram
/// which is then saved to file (appending .hist to the filename),
/// also printing histogram statistics.
///
/// The resulting histogram file can be easily be visualized using gnuplot:
/// plot "filename.hist" with boxes
///

#include "theoretica.h"
using namespace th;

#include <fstream>
#include <iostream>


int main(int argc, const char** argv) {

    if (argc < 2) {
        std::cout << "Usage: histogram <filename>" << std::endl;
        return 1;
    }

    std::string filename = argv[1];
    std::ifstream infile (filename);

    if (!infile.is_open()) {
        std::cout << "Unable to open input file." << std::endl;
        return 2;
    }

    // Construct histogram while reading from file
    std::vector<real> data;
    std::string line;

    std::cout << "Reading data points from file..." << std::endl;

    while(std::getline(infile, line)) {

        if(line == "")
            continue;

        data.push_back(std::stod(line));
    }

    std::ofstream outfile (filename + ".hist");

    if (!outfile.is_open()) {
        std::cout << "Unable to open output file." << std::endl;
        return 3; 
    }

    std::cout << "Constructing histogram from data..." << std::endl;

    unsigned int bins = int(sqrt(data.size()));
    if (argc > 2)
        bins = std::stoi(argv[2]);

    histogram h (data, bins);
    outfile << h;

    std::cout << "Wrote histogram to file" << std::endl;

    // Print to standard output some useful statistics
    std::cout << "Statistics:" << std::endl;
    std::cout << "N = " << h.number() << std::endl; 
    std::cout << "Mean: " << stats::mean(h) << std::endl; 
    std::cout << "Variance: " << stats::variance(h) << std::endl; 
    std::cout << "Standard Deviation: " << stats::stdev(h) << std::endl;
}
