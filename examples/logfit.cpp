
///
/// @file logfit.cpp Logarithmic fit from file data
///

#include "theoretica.h"
using namespace th;

#include <iostream>
#include <fstream>
#include <sstream>


int main(int argc, char const *argv[]) {
    
    if(argc < 2) {
        std::cout << "Usage: logfit <filename>" << std::endl;
        return 1;
    }

    std::string filename = std::string(argv[1]);
    std::ifstream file(filename);

    if(!file.is_open()) {
        std::cout << "Unable to open file" << std::endl;
        return 2;
    }

    std::string line;
    vec_buff X;
    vec_buff Y;

    while(std::getline(file, line)) {
        
        real x, y;
        std::stringstream s(line);

        s >> x;
        s >> y;

        X.push_back(th::ln(x));
        Y.push_back(th::ln(y));
    }

    real intercept = lst_sqrs_lin_intercept(X, Y);
    real slope = lst_sqrs_lin_slope(X, Y);

    std::cout << "Model: y = A + Bx" << std::endl;
    std::cout << "A = " << intercept << std::endl;
    std::cout << "B = " << slope << std::endl;
    std::cout << "Err. = " << lst_sqrs_lin_error(X, Y, intercept, slope) << std::endl;

    return 0;
}
