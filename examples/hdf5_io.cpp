

#ifdef THEORETICA_HAS_HDF5

#include "theoretica.h"
#include "io/hdf5.h"

using namespace theoretica;


int main() {

	try {

		io::println("Creating HDF5 file...");
		io::hdf5_file file ("file.h5", true);

		io::println("\nFile contents:");
		io::println(file);

		io::println("Creating group...");
		file.create_group("/group");

		io::println("Writing vector to file...");
		vec<real> v (10000);
		PRNG g (151315991);
		for (real& x : v)
			x = rand_gaussian(0, 1, g);
		file.write_vec("/group/vec", v);

		mat<real> A (1000, 1000);
		for (real& x : A)
			x = rand_gaussian(0, 1, g);
		file.write_mat("/group/mat", A);

		file.write_attribute("/group/vec", "author", std::string("Albert Einstein"));
		file.write_attribute("/group/mat", "desc", std::string("A matrix of 1000x1000 random Gaussian values"));

		file.refresh();

		io::println("\nAttributes:");
		io::println("Author =", file.read_attribute<std::string>("/group/vec", "author"));
		io::println("Description =", file.read_attribute<std::string>("/group/mat", "desc"));

		io::println("\nFile contents:");
		io::println(file);

	} catch (const std::exception& e) {
		io::println("Error creating HDF5 file: " + std::string(e.what()));
		return 1;
	}
	
	return 0;
}

#else

// Stub for when HDF5 support is not enabled
#include <iostream>
int main() {
	std::cout << "HDF5 support is not enabled. To enable it, set HDF5_INCLUDE and HDF5_LIB variables in the Makefile." << std::endl;
	return 0;
}

#endif