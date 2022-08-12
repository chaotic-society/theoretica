
///
/// @file utility.h Utility functions
/// This header file is **not** automatically included by theoretica.h
///

#ifndef THEORETICA_UTILITY_H
#define THEORETICA_UTILITY_H

#include <iostream>
#include <string>
#include <algorithm>

#include "./core/vec_buff.h"


// Log a variable to stdout
#define LOG(x) (std::cout << #x << " = " << x << std::endl)


namespace theoretica {


	/// Insert a `vec_buff` data set from standard input
	inline void insert_data(vec_buff& data, std::string terminator = "") {

		std::string line;
		real value;

		while(true) {
			std::getline(std::cin, line);

			if(line == terminator)
				break;

			if(line == "")
				continue;

			std::replace(line.begin(), line.end(), ',', '.');

			try {
				value = std::stod(line);
			} catch(...) {
				std::cout << "Input conversion error" << std::endl;
				value = nan();
			}

			if(!is_nan(value))
				data.emplace_back(value);
		}
	}

}

#endif
