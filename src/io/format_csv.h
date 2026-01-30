
///
/// @file format_csv.h CSV file format support.
///

#ifndef THEORETICA_FORMAT_CSV_H
#define THEORETICA_FORMAT_CSV_H

#include <fstream>

#include "../algebra/vec.h"
#include "../algebra/mat.h"


namespace theoretica {
namespace io {

	
	/// Write a vector to file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param v The vector to write
	/// @param header An optional header to add for the vector column
	template<typename Type, unsigned int N>
	inline void write_csv(const std::string& filename, const vec<Type, N>& v, std::string header = "") {

		std::ofstream file (filename);

		if (!file.is_open()) {
			
			// TODO: throw another exception ?

			TH_MATH_ERROR("io::write_csv", false, MathError::ImpossibleOperation);
			return;
		}

		if (header != "")
			file << "\"" << header << "\"" << std::endl;

		for (size_t i = 0; i < v.size(); ++i)
			file << v[i] << std::endl;
	}


	/// Read a vector from a file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param v A reference to the vector to overwrite
	template<typename Type, unsigned int N>
	inline void read_csv(const std::string& filename, vec<Type, N>& v) {

	}


	/// Write a matrix to file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param A The matrix to write
	template<typename Type, unsigned int N, unsigned int M>
	inline void write_csv(const std::string& filename, const mat<Type, N, M>& A, const std::string& separator = ", ") {

		std::ofstream file (filename);

		if (!file.is_open()) {
			
			// TODO: throw another exception ?

			TH_MATH_ERROR("io::write_csv", false, MathError::ImpossibleOperation);
			return;
		}

		//if (header != "") file << "\"" << header << "\"" << std::endl;

		for (size_t i = 0; i < A.rows(); i++) {
			for (size_t j = 0; j < A.cols(); j++) {
				
				file << A(i, j);

				if (j != A.cols() - 1)
					file << separator;
				else
					file << std::endl;
			}
		}
	} 


}}

#endif
