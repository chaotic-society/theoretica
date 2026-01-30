
///
/// @file data_table.h Data table structure for holding labeled columns of data.
///

#ifndef THEORETICA_DATA_TABLE_H
#define THEORETICA_DATA_TABLE_H

#include <map>
#include <string>

#include "../core/constants.h"
#include "../algebra/vec.h"


namespace theoretica {

	/// @class data_table
	/// A data structure holding a table of numerical values,
	/// such as experimental datasets, where each column has a name.
	class data_table {
	private:

		// TO-DO: Change to std::variant over vectors for multiple data types ?

		/// Hashmap holding the columns by column name
		std::map<std::string, vec<real>> table;

		/// Ordered list of column names
		std::vector<std::string> indicies;

	public:

		/// Default constructor
		data_table() {}


		/// Construct from a list of the names
		/// of the column datasets.
		data_table(std::initializer_list<std::string> l) {

			indicies.reserve(l.size());
			size_t i = 0;

			for (const auto& name : l) {

				table.insert(std::make_pair(name, vec<real>()));
				
				indicies[i] = name;
				i++;
			}
		}


		/// Get a reference to a column by column index.
		inline vec<real>& operator[](unsigned int i) {
			return table[indicies[i]];
		}


		/// Get a reference to a column by column name.
		inline vec<real>& operator[](const std::string& name) {
			return table[name];
		}


		/// Insert a new column in the last place.
		inline void insert(std::string name, const vec<real>& v) {

			table[name] = v;

			// Column doesn't exist yet.
			if(table.find(name) == table.end())
				indicies.emplace_back(name);
		}
		
	};

}


#endif
