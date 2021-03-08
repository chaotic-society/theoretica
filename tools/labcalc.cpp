#define UROBORO_INCLUDE_ALL
#include "./statistics.h"
#include "./utility.h"
#include <iostream>
#include <string>

using namespace uroboro;

int main(int argc, char const *argv[]) {

	std::cout << "Insert each value and press Enter, write END to stop insertion" << std::endl;

	vec_buff X;
	std::string line;
	real value;

	while(true) {

		std::getline(std::cin, line);

		if(line == "END")
			break;

		try {
			value = std::stod(line);
		} catch(...) {
			std::cout << "Input conversion error" << std::endl;
		}

		X.emplace_back(value);

	}

	print_sample_stats(X);

	std::cout << "Press Enter to exit..." << std::endl;
	std::cin.get();

	return 0;
}
