
///
/// @file data_table.h Data table structure for holding labeled columns of data.
///

#ifndef THEORETICA_DATA_TABLE_H
#define THEORETICA_DATA_TABLE_H

#include <map>
#include <vector>
#include <unordered_map>
#include <string>
#include <iomanip>

#include "../core/constants.h"
#include "../algebra/vec.h"
#include "../algebra/mat.h"


namespace theoretica {


	/// @class data_table
	/// A data structure for holding labeled columns of data,
	/// where each column is a vector of real numbers.
	class data_table {
	private:

		// Vector of data columns.
		std::vector<vec<real>> columns {};

		// Ordered list of column names
		std::vector<std::string> column_names {};

		// Map of column names to their indices in the data table
		std::unordered_map<std::string, size_t> indices {};

	public:

		/// Default constructor
		data_table() {}


		/// Construct from map of column vectors.
		/// The order of columns is determined by the order of the map iteration.
		///
		/// @param table A map where keys are column names and values are column vectors.
		data_table(const std::map<std::string, vec<real>>& table) {

			for (const auto& pair : table) {
				column_names.emplace_back(pair.first);
				columns.emplace_back(pair.second);
				indices[pair.first] = columns.size() - 1;
			}
		}
		

		/// Construct a data table with preallocated size
		///
		/// @param num_rows The number of rows in the data table
		/// @param column_names The names of the columns in the data table
		data_table(size_t num_rows, const std::vector<std::string>& column_names) {

			this->column_names = column_names;
			this->columns.resize(column_names.size(), vec<real>(num_rows));

			for (size_t i = 0; i < column_names.size(); ++i)
				indices[column_names[i]] = i;
		}
		

		/// Copy constructor from another data table.
		data_table(const data_table& other) {
			columns = other.columns;
			column_names = other.column_names;
			indices = other.indices;
		}

		/// Move constructor from another data table, leaving the other empty.
		data_table(data_table&& other) noexcept {

			columns = std::move(other.columns);
			column_names = std::move(other.column_names);
			indices = std::move(other.indices);
			other.columns.clear();
			other.column_names.clear();
			other.indices.clear();
		}


		/// Get the (maximum) number of rows in the data table.
		/// This operation is O(K) in the number of columns K.
		///
		/// @return The size of the first column or 0 if the table is empty.
		inline size_t rows() const {

			if (columns.empty())
				return 0;

			unsigned int max_rows = 0;
			for (const auto& col : columns)
				if (col.size() > max_rows)
					max_rows = col.size();

			return max_rows;
		}


		/// Get the number of columns in the data table.
		/// @return The number of columns in the data table.
		inline size_t cols() const {
			return columns.size();
		}


		/// Get the total number of elements in the data table, i.e. N. rows * N. columns.
		/// @return The total number of elements in the data table.
		inline size_t size() const {

			size_t sum = 0;
			for (const auto& col : columns)
				sum += col.size();

			return sum;
		}


		/// Check whether the data table is empty (i.e. has no columns or all columns are empty).
		/// @return True if the data table is empty, false otherwise.
		inline bool empty() const {
			return size() == 0;
		}


		/// Remove all columns from the data table, leaving it empty.
		inline void clear() {
			columns.clear();
			column_names.clear();
			indices.clear();
		}


		/// Get a list of the column names in the data table.
		/// @return A vector of column names in the data table.
		inline std::vector<std::string> header() const {
			return column_names;
		}


		/// Check whether the data table has a column with the given name.
		/// @param name The name of the column to check for.
		/// @return Whether the table contains a column with the given name.
		inline bool has_column(const std::string& name) const {
			return indices.find(name) != indices.end();
		}


		/// Get the columns of the data table as a simple vector of column vectors,
		/// without the column names or indices, by reference.
		///
		/// @return A vector of column vectors in the data table.
		inline std::vector<vec<real>>& data() {
			return columns;
		}


		/// Get the columns of the data table as a simple vector of column vectors, without the column names or indices.
		///
		/// @return A vector of column vectors in the data table.
		inline const std::vector<vec<real>>& data() const {
			return columns;
		}


		/// Access a column by name, returning a reference to the column vector.
		/// Throws an out_of_range exception if the column name is not found.
		///
		/// @param name The name of the column to access.
		/// @return A reference to the column vector with the given name.
		inline vec<real>& operator[](const std::string& name) {

			auto it = indices.find(name);
			if (it == indices.end()) {
				throw std::out_of_range("data_table::operator[]: column '" + name + "' not found");
			}

			return columns[it->second];
		}


		/// Access a column by name, returning a const reference to the column vector.
		/// Throws an out_of_range exception if the column name is not found.
		/// @param name The name of the column to access.
		/// @return A const reference to the column vector with the given name.
		inline const vec<real>& operator[](const std::string& name) const {

			auto it = indices.find(name);
			if (it == indices.end()) {
				throw std::out_of_range("data_table::operator[]: column '" + name + "' not found");
			}

			return columns[it->second];
		}


		/// Access a column by name, returning a reference to the column vector.
		/// Throws an out_of_range exception if the column name is not found.
		/// @param name The name of the column to access.
		/// @return A reference to the column vector with the given name.
		inline vec<real>& at(const std::string& name) {

			auto it = indices.find(name);
			if (it == indices.end()) {
				throw std::out_of_range("data_table::at: column '" + name + "' not found");
			}

			return columns[it->second];
		}


		/// Access a column by index, returning a reference to the column vector.
		/// Does not perform bounds checking, so the behavior is undefined if the index is out of range,
		/// for safer access consider using the at() method instead, which performs bounds checking.
		/// @param idx The index of the column to access.
		/// @return A reference to the column vector at the given index.
		inline vec<real>& operator[](size_t idx) {
			return columns[idx];
		}


		/// Access a column by index, returning a const reference to the column vector.
		/// Does not perform bounds checking, so the behavior is undefined if the index is out of range,
		/// for safer access consider using the at() method instead, which performs bounds checking.
		/// @param idx The index of the column to access.
		/// @return A const reference to the column vector at the given index.
		inline const vec<real>& column(size_t idx) const {
			return columns[idx];
		}


		/// Select a subset of columns from the table, returning a new table containing only the selected columns.
		/// If a column name in the selection list is not found in the data table, it is ignored.
		///
		/// @note The current implementation immediately returns a new data table,
		/// but in the future it may be optimized to return a view of the original data
		/// table instead, to avoid unnecessary copying.
		///
		/// @param cols A vector of column names to select from the data table.
		/// @return A new data table containing only the selected columns.
		inline data_table select(const std::vector<std::string>& cols) const {

			data_table result;

			for (const auto& col_name : cols) {

				auto it = indices.find(col_name);
				if (it != indices.end()) {
					result.column_names.emplace_back(col_name);
					result.columns.emplace_back(columns[it->second]);
					result.indices[col_name] = result.columns.size() - 1;
				}
			}

			return result;
		}


		/// Access an element by row index and column name, returning a reference to the value.
		/// Throws an out_of_range exception if the column name is not found or the row index is out of range.
		///
		/// @param row The index of the row to access.
		/// @param col The name of the column to access.
		/// @return A reference to the value at the given row and column.
		inline real& at(const std::string& col, size_t row) {

			auto it = indices.find(col);
			if (it == indices.end()) {
				throw std::out_of_range("data_table::at: column '" + col + "' not found");
			}

			vec<real>& column_vec = columns[it->second];
			if (row >= column_vec.size()) {
				throw std::out_of_range("data_table::at: row index " + std::to_string(row) + " out of range");
			}

			return column_vec[row];
		}


		/// Access an element by row index and column name, returning a const reference to the value.
		/// Throws an out_of_range exception if the column name is not found or the row index is out of range.
		///
		/// @param row The index of the row to access.
		/// @param col The name of the column to access.
		/// @return A const reference to the value at the given row and column.
		inline const real& at(const std::string& col, size_t row) const {
			
			auto it = indices.find(col);
			if (it == indices.end()) {
				throw std::out_of_range("data_table::at: column '" + col + "' not found");
			}

			const vec<real>& column_vec = columns[it->second];
			if (row >= column_vec.size()) {
				throw std::out_of_range("data_table::at: row index " + std::to_string(row) + " out of range");
			}

			return column_vec[row];
		}


		/// Get an entire row as a map of column names to values.
		/// Missing values for columns that do not have enough rows are filled with NaN.
		///
		/// @param idx The index of the row to access.
		/// @return A map of column names to values for the given row.
		inline std::unordered_map<std::string, real> row(size_t idx) const {

			std::unordered_map<std::string, real> result;

			for (size_t i = 0; i < columns.size(); ++i) {
				const real value = idx < columns[i].size() ? columns[i][idx] : nan();
				result[column_names[i]] = value;
			}

			return result;
		}

		
		/// Get an entire row as a vector of values, in the same order as the columns in the data table.
		/// Missing values for columns that do not have enough rows are filled with NaN.
		///
		/// @param idx The index of the row to access.
		/// @return A map of column names to values for the given row.
		inline vec<real> row_vec(size_t idx) const {

			vec<real> result(columns.size());

			for (size_t i = 0; i < columns.size(); ++i)
				result[i] = idx < columns[i].size() ? columns[i][idx] : nan();

			return result;
		}


		/// Get the first n rows of the data table as a new data table.
		/// If n is greater than the number of rows in the table, the entire table is returned.
		///
		/// @param n The number of rows to include in the head of the data table.
		/// @return A new data table containing the first n rows of the original table.
		inline data_table head(size_t n = 5) const {

			data_table result;

			for (size_t i = 0; i < columns.size(); ++i) {

				vec<real> column_head (min(n, columns[i].size()));

				for (size_t j = 0; j < column_head.size(); ++j)
					column_head[j] = columns[i][j];

				result.column_names.emplace_back(column_names[i]);
				result.columns.emplace_back(column_head);
				result.indices[column_names[i]] = result.columns.size() - 1;
			}

			return result;
		}


		/// Get the last n rows of the data table as a new data table.
		/// If n is greater than the number of rows in the table, the entire table is returned.
		/// @param n The number of rows to include in the tail of the data table.
		/// @return A new data table containing the last n rows of the original table.
		inline data_table tail(size_t n = 5) const {
			
			data_table result;

			for (size_t i = 0; i < columns.size(); ++i) {

				vec<real> column_tail(min(n, columns[i].size()));

				for (size_t j = 0; j < column_tail.size(); ++j)
					column_tail[j] = columns[i][columns[i].size() - column_tail.size() + j];

				result.column_names.emplace_back(column_names[i]);
				result.columns.emplace_back(column_tail);
				result.indices[column_names[i]] = result.columns.size() - 1;
			}

			return result;
		}


		/// Insert a new column into the data table with the given name and data.
		/// If a column with the same name already exists, it is overwritten.
		///
		/// @param name The name of the column to add.
		/// @param data The data representing the new column.
		inline data_table& insert(const std::string& name, const vec<real>& data) {

			auto it = indices.find(name);
			if (it != indices.end()) {
				columns[it->second] = data;
			} else {
				column_names.emplace_back(name);
				columns.emplace_back(data);
				indices[name] = columns.size() - 1;
			}

			return *this;
		}


		/// Insert a new column into the data table with the given name, number of rows, and constant value.
		/// If a column with the same name already exists, it is overwritten.
		///
		/// @param name The name of the column to add.
		/// @param num_rows The number of rows in the new column.
		/// @param value The value to fill the new column with.
		inline data_table& insert(const std::string& name, size_t num_rows, real value = 0.0) {

			vec<real> data (num_rows, value);
			insert(name, data);

			return *this;
		}


		/// Drop a column from the data table by name. If the column name is not found, the table is unchanged.
		/// This operation is O(K) in the number of columns K, since it requires updating the indices of all subsequent columns.
		///
		/// @param name The name of the column to drop.
		inline data_table& drop_column(const std::string& name) {

			auto it = indices.find(name);
			if (it != indices.end()) {

				size_t idx = it->second;
				columns.erase(columns.begin() + idx);
				column_names.erase(column_names.begin() + idx);
				indices.erase(it);

				// Update indices of remaining columns
				for (size_t i = idx; i < column_names.size(); ++i)
					indices[column_names[i]] = i;
			}

			return *this;
		}


		/// Drop multiple columns from the data table by name. Column names that are not found are ignored.
		/// This operation is O(K) with respect to the number of columns K, since it requires
		/// updating the indices of all subsequent columns for each dropped column.
		///
		/// @param names The names of the columns to drop.
		inline data_table& drop_columns(const std::vector<std::string>& names) {

			for (const auto& name : names)
				drop_column(name);

			return *this;
		}


		/// Rename a column in the data table from old_name to new_name. If old_name is not found, the table is unchanged,
		/// while if new_name already exists, it is overwritten.
		///
		/// @param old_name The old name of the column to rename
		/// @param new_name The new name of the column
		inline data_table& rename(const std::string& old_name, const std::string& new_name) {

			auto it = indices.find(old_name);
			if (it != indices.end()) {

				size_t i = it->second;
				column_names[i] = new_name;
				indices.erase(it);
				indices[new_name] = i;
			}

			return *this;
		}


		/// Convert the data table to a matrix, where each column of the matrix
		/// corresponds to a column in the data table, and each row corresponds to
		/// a row in the data table. Missing values are filled with NaN.
		///
		/// @return A matrix representation of the data table.
		inline mat<real> to_matrix() const {
			
			mat<real> result (columns.empty() ? 0 : columns[0].size(), columns.size());

			for (size_t i = 0; i < columns.size(); ++i) {

				for (size_t j = 0; j < result.rows(); ++j)
					result(j, i) = j < columns[i].size() ? columns[i][j] : nan();
			}

			return result;
		}


		/// Convert a matrix to a data table, where each column of the matrix corresponds to a column in the data table,
		/// and each row corresponds to a row in the data table. The column names are taken from the provided vector,
		/// and if there are fewer column names than columns in the matrix, the remaining columns
		/// are named "col_i" where i is the column index. If there are more column names than columns in the matrix,
		/// the extra column names are ignored.
		///
		/// @param m The matrix to read elements from.
		/// @param col_names Names of the columns in the resulting data table.
		inline data_table& from_matrix(const mat<real>& m, const std::vector<std::string>& col_names) {

			columns.clear();
			column_names.clear();
			indices.clear();

			for (size_t i = 0; i < m.cols(); ++i) {

				vec<real> column(m.rows());

				for (size_t j = 0; j < m.rows(); ++j)
					column[j] = m(j, i);

				std::string col_name = i < col_names.size() ? col_names[i] : "col_" + std::to_string(i);
				column_names.emplace_back(col_name);
				columns.emplace_back(column);
				indices[col_name] = columns.size() - 1;
			}

			return *this;
		}


#ifndef THEORETICA_NO_PRINT

		/// Convert the data table to string representation
		inline std::string to_string(unsigned int max_rows = 8, unsigned int precision = 6, unsigned int max_width = 12) const {

			std::stringstream res;

			// Print column names
			for (size_t i = 0; i < column_names.size(); ++i) {

				if (column_names[i].size() > max_width)
					res << std::setw(max_width + 2) << column_names[i].substr(0, max_width - 3) + "...";
				else
					res << std::setw(max_width + 2) << column_names[i] << "\t";
			}

			res << "\n";

			// Print data rows
			size_t num_rows = columns.empty() ? 0 : columns[0].size();
			size_t rows_to_print = min(max_rows, num_rows);

			for (size_t i = 0; i < rows_to_print; ++i) {

				for (size_t j = 0; j < columns.size(); ++j) {

					if (i < columns[j].size()) {
						res << std::setw(max_width + 2) << std::setprecision(precision) << columns[j][i] << "\t";
					} else {
						res << std::setw(max_width + 2) << "\t";
					}
				}

				res << "\n";
			}

			if (num_rows > max_rows) {
				res << "... " << (num_rows - max_rows) << " more rows\n";
			}

			return res.str();
		}


		/// Convert the data table to string representation
		inline operator std::string() {
			return to_string();
		}


		/// Stream the data table in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(std::ostream& out, const data_table& obj) {
			return out << obj.to_string();
		}

#endif
		
	};

}


#endif
