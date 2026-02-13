
///
/// @file format_csv.h CSV file format support.
///

#ifndef THEORETICA_FORMAT_CSV_H
#define THEORETICA_FORMAT_CSV_H

#include <fstream>
#include <iomanip>

#include "../algebra/vec.h"
#include "../algebra/mat.h"
#include "./strings.h"


namespace theoretica {
namespace io {

	
	/// Write a vector to file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param v The vector to write
	/// @param header An optional header to add for the vector column
	template<typename Type, unsigned int N>
	inline void write_csv(
		const std::string& filename, const vec<Type, N>& v,
		const std::string& header = "", unsigned int precision = 8) {

		std::ofstream file (filename);

		if (!file.is_open()) {
			
			// TODO: throw another exception ?
			TH_MATH_ERROR("io::write_csv", false, MathError::ImpossibleOperation);
			return;
		}

		if (header != "")
			file << "\"" << header << "\"" << std::endl;

		for (size_t i = 0; i < v.size(); ++i)
			file << std::setprecision(precision) << v[i] << std::endl;
	}


	/// Read a vector from a file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param v A reference to the vector to overwrite
	template<unsigned int N>
	inline void read_csv(const std::string& filename, vec<real, N>& v) {

		std::ifstream file (filename);
		std::string line;

		if (!file.is_open()) {
			
			// TODO: throw another exception ?
			TH_MATH_ERROR("io::read_csv", false, MathError::ImpossibleOperation);
			return;
		}

		// Check for header
		std::getline(file, line);
		line = io::trim(line);

		// Resulting column vector
		std::vector<real> col;
		
		if (io::is_number(line)) {
			try {
				real first = std::stod(line);
				col.emplace_back(first);
			} catch (const std::invalid_argument& e) {}
		}

		while (std::getline(file, line)) {

			line = io::trim(line);

			try {
				real val = std::stod(line);
				col.emplace_back(val);
			} catch (const std::invalid_argument& e) {
				col.emplace_back(nan());
			}
		}

		// Handle mismatched sizes with empty values (NaN)
		if (v.size() > col.size()) {

			for (size_t i = 0; i < col.size(); i++)
				v[i] = col[i];

			for (size_t i = col.size(); i < v.size(); i++)
				v[i] = nan();
		} else {
			v = col;
		}
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
