
///
/// @file format_csv.h CSV file format support.
///

#ifndef THEORETICA_FORMAT_CSV_H
#define THEORETICA_FORMAT_CSV_H

#include <fstream>
#include <iomanip>
#include <algorithm>

#include "../algebra/vec.h"
#include "../algebra/mat.h"
#include "./data_table.h"
#include "./strings.h"


namespace theoretica {
namespace io {


	/// Parse a CSV line handling quoted fields.
	/// Supports fields enclosed in double quotes, with delimiters (commas) inside quoted fields.
	/// Does not support escaped quotes or multiline fields.
	///
	/// @param line The CSV line to parse
	/// @param delimiter The field delimiter (default is comma)
	/// @return A vector of parsed fields
	inline std::vector<std::string> parse_csv(const std::string& line, char delimiter = ',') {

		std::vector<std::string> fields;
		std::string field;
		bool quoted = false;

		for (size_t i = 0; i < line.length(); ++i) {

			char c = line[i];

			if (c == '"') {
				quoted = !quoted;
			} else if (c == delimiter && !quoted) {
				fields.emplace_back(field);
				field.clear();
			} else if(!std::isspace(c) || quoted) {
				field += c;
			}
		}

		fields.emplace_back(field);
		return fields;
	}


	/// Given a string entry, sanitize it for printing to a CSV file.
	/// If it contains commas or whitespace, quotes are added before and after the string.
	///
	/// @param str The string to check for relevant characters.
	/// @return The sanitized string, ready for printing to CSV.
	inline std::string quote_csv(const std::string& str) {

		bool has_whitespace = false;
		for (char c : str) {

			if (std::isspace(c)) {
				has_whitespace = true;
				break;
			}
		}

		if (str.find(',') != std::string::npos || has_whitespace)
			return "\"" + str + "\"";
		else
			return str;
	}

	
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
	/// If present, the header is ignored. If the file has fewer elements than the vector,
	/// the remaining elements are filled with NaN.
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
		if (!std::getline(file, line))
			return;

		line = io::unquote(io::trim(line));

		// Resulting column vector
		std::vector<real> col;

		if (io::is_number(line)) {
			try {
				std::replace(line.begin(), line.end(), ',', '.');
				real first = std::stod(line);
				col.emplace_back(first);
			} catch (const std::invalid_argument& e) {}
		}

		// All remaining lines are data
		while (std::getline(file, line)) {

			line = io::unquote(io::trim(line));
			std::replace(line.begin(), line.end(), ',', '.');

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
	inline void write_csv(
		const std::string& filename, const mat<Type, N, M>& A,
		const std::string& delimiter = ", ", unsigned int precision = 8) {

		std::ofstream file (filename);

		if (!file.is_open()) {
			
			// TODO: throw another exception ?
			TH_MATH_ERROR("io::write_csv", false, MathError::ImpossibleOperation);
			return;
		}

		for (size_t i = 0; i < A.rows(); i++) {
			for (size_t j = 0; j < A.cols(); j++) {
				
				file << std::setprecision(precision) << A(i, j);

				if (j != A.cols() - 1)
					file << delimiter;
				else
					file << std::endl;
			}
		}
	}


	/// Read a matrix from a file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param A A reference to the matrix to overwrite
	/// If present, the header is ignored. If the file has fewer rows or columns than the matrix,
	/// the remaining elements are filled with NaN.
	template<unsigned int N, unsigned int K>
	inline void read_csv(const std::string& filename, mat<real, N, K>& A) {

		std::ifstream file (filename);
		std::string line;

		if (!file.is_open()) {
			
			// TODO: throw another exception ?
			TH_MATH_ERROR("io::read_csv", false, MathError::ImpossibleOperation);
			return;
		}

		std::vector<std::vector<real>> rows;

		// Read first line to check for header
		if (!std::getline(file, line))
			return;

		std::vector<std::string> first_row = parse_csv(line);

		// Check if first line is a header
		bool has_header = false;

		for (const auto& cell : first_row) {

			if (!io::is_number(cell)) {
				has_header = true;
				break;
			}
		}

		// If first line is not a header, process it as data
		if (!has_header) {

			std::vector<real> row;

			for (auto cell : first_row) {

				std::replace(cell.begin(), cell.end(), ',', '.');

				try {
					row.emplace_back(std::stod(cell));
				} catch (const std::invalid_argument& e) {
					row.emplace_back(nan());
				}
			}

			if (!row.empty())
				rows.emplace_back(row);
		}

		// Read remaining lines
		while (std::getline(file, line)) {

			// Skip empty lines
			if (line.empty())
				continue;

			std::vector<std::string> cells = parse_csv(line);
			std::vector<real> row;

			for (auto cell : cells) {

				std::replace(cell.begin(), cell.end(), ',', '.');

				try {
					row.emplace_back(std::stod(cell));
				} catch (const std::invalid_argument& e) {
					row.emplace_back(nan());
				}
			}

			if (!row.empty()) {
				rows.emplace_back(row);
			}
		}

		A.resize(rows[0].size(), rows.size());

		if (A.rows() < rows[0].size() || A.cols() < rows.size()) {
			TH_MATH_ERROR("io::read_csv", false, MathError::ImpossibleOperation);
			algebra::mat_error(A);
			return;
		}

		// Fill matrix with parsed data
		for (size_t i = 0; i < min(rows.size(), A.rows()); ++i) {

			for (size_t j = 0; j < min(rows[i].size(), A.cols()); ++j)
				A(i, j) = rows[i][j];

			// Pad remaining columns with NaN
			for (size_t j = rows[i].size(); j < A.cols(); ++j)
				A(i, j) = nan();
		}

		// Pad remaining rows with NaN
		for (size_t i = rows.size(); i < A.rows(); ++i)
			for (size_t j = 0; j < A.cols(); ++j)
				A(i, j) = nan();
	}


	/// Write a data_table to file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param table The data_table to write
	inline void write_csv(
		const std::string& filename, const data_table& table,
		const std::string& delimiter = ", ", unsigned int precision = 8) {

		std::ofstream file (filename);

		if (!file.is_open()) {
			TH_MATH_ERROR("io::write_csv", false, MathError::ImpossibleOperation);
			return;
		}

		// Write header
		bool first = true;
		for (const std::string& name : table.header()) {

			if (!first)
				file << delimiter;

			file << quote_csv(name);
			first = false;
		}
		file << std::endl;

		// Write data rows
		size_t max_rows = table.rows();
		for (size_t i = 0; i < max_rows; ++i) {

			first = true;
			for (const auto& col : table.data()) {

				if (!first)
					file << delimiter;

				if (i < col.size())
					file << std::setprecision(precision) << col[i];
				else
					file << nan();

				first = false;
			}
			file << std::endl;
		}

	}


	/// Read a data_table from a file in the CSV format.
	///
	/// @param filename The name of the file
	/// @param table A reference to the data_table to overwrite
	/// If present, the header is read as column names. If the file has fewer rows than the data_table,
	/// the remaining rows are filled with NaN.
	inline void read_csv(const std::string& filename, data_table& table) {

		std::ifstream file (filename);
		std::string line;

		if (!file.is_open()) {
			TH_MATH_ERROR("io::read_csv", false, MathError::ImpossibleOperation);
			return;
		}

		// Read header
		if (!std::getline(file, line))
			return;

		std::vector<std::string> first_row = parse_csv(line);
		
		if (first_row.empty())
			return;
		
		std::vector<std::string> column_names;
		size_t num_cols = first_row.size();
		std::vector<vec<real>> columns (num_cols);

		// Check if first line is a header
		bool has_header = false;

		for (const auto& cell : first_row) {

			if (!io::is_number(cell)) {
				has_header = true;
				break;
			}
		}

		// If first line is not a header, process it as data
		if (!has_header) {

			for (size_t j = 0; j < first_row.size(); ++j) {

				std::string cell = first_row[j];
				std::replace(cell.begin(), cell.end(), ',', '.');

				try {
					columns[j].append(std::stod(cell));
				} catch (const std::invalid_argument& e) {
					columns[j].append(nan());
				}
			}
		
			// Generate default column names
			for (size_t j = 0; j < num_cols; ++j) {
				column_names.emplace_back("col" + std::to_string(j));
			}
		
		} else {
		
			for (const auto& name : first_row)
				column_names.emplace_back(io::unquote(io::trim(name)));
		}
		first_row.clear();

		// Read data rows
		while (std::getline(file, line)) {

			if (line.empty())
				continue;

			std::vector<std::string> cells = parse_csv(line);

			for (size_t j = 0; j < num_cols; ++j) {

				if (j < cells.size()) {
					std::string cell = cells[j];
					std::replace(cell.begin(), cell.end(), ',', '.');

					try {
						real val = std::stod(cell);
						columns[j].append(val);
					} catch (const std::invalid_argument& e) {
						columns[j].append(nan());
					}
				} else {
					columns[j].append(nan());
				}
			}
		}
		
		for (size_t j = 0; j < num_cols; ++j)
			table.insert(column_names[j], columns[j]);
	}


	/// Read a generic data structure from a file in the CSV format,
	/// specifying the target type. Supported types include vec<> and mat<>.
	///
	/// @param filename The name of the file
	/// @tparam Type The type of the data structure to read
	/// @return The data structure read from the file
	/// The template parameter Type should specify
	/// a valid supported type, e.g. mat<real>, mat3, vec<real> and so on:
	/// auto A = read_csv<mat3>("myfile.csv");
	template<typename Type>
	inline Type read_csv(const std::string& filename) {
		Type A;
		read_csv(filename, A);
		return A;
	}

}}

#endif
