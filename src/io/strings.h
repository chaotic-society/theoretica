
#ifndef THEORETICA_STRINGS_H
#define THEORETICA_STRINGS_H

#include <string>


namespace theoretica {
namespace io {


	/// Check if a given string could be correctly interpreted as a number.
	///
	/// @param str The string to check
	/// @return Whether the string could be interpreted as a number
	inline bool is_number(const std::string& str) {

		if (str.empty())
			return false;

		const std::string allowed = "1234567890.,Ee+-";

		for (char c : str)
			if (allowed.find(c) == std::string::npos)
				return false;

		return true;
	}


	/// Remove all leading and trailing whitespace from a string,
	/// returning the resulting string.
	///
	/// @param str The input string
	/// @return The trimmed string
	inline std::string trim(const std::string& str) {
		
		size_t start = 0;
		while (start < str.length() && std::isspace(str[start]))
			++start;
		
		size_t end = str.length();
		while (end > start && std::isspace(str[end - 1]))
			--end;
		
		return str.substr(start, end - start);
	}


	/// Remove leading and trailing double quotes from a string, if present.
	///
	/// @param str The input string
	/// @return The string with quotes removed
	inline std::string unquote(const std::string& str) {

		if (str.length() >= 2 && str.front() == '"' && str.back() == '"')
			return str.substr(1, str.length() - 2);

		return str;
	}

}}


#endif
