
///
/// @file utility.h Optional header with utilities for input and output.
///

#ifndef THEORETICA_UTILITY_H
#define THEORETICA_UTILITY_H

#include <iostream>
#include <string>
#include <algorithm>
#include "./algebra/vec.h"


namespace theoretica {


	/// Print the given argument to standard output.
	template<typename Type>
	inline void print(const Type& curr) {
		std::cout << curr;
	}

	/// Print the given arguments to standard output
	/// separated by a space.
	template<typename Type, typename ...Args>
	inline void print(const Type& curr, Args... args) {
		std::cout << curr << " ";
		print(args...);
	}

	/// Print the given argument to a stream.
	template<typename Type>
	inline void fprint(std::ostream& out, const Type& last) {
		out << last;
	}

	/// Print the given arguments to standard output
	/// separated by a space.
	template<typename Type, typename ...Args>
	inline void fprint(std::ostream& out, const Type& curr, Args... args) {
		out << curr << " ";
		fprint(out, args...);
	}


	/// Print a newline to standard output.
	inline void println() {
		std::cout << "\n";
	}

	/// Print the given argument to standard output
	/// followed by a newline.
	template<typename Type>
	inline void println(const Type& curr) {
		std::cout << curr << "\n";
	}

	/// Print the given arguments to standard output
	/// separated by a space and followed by a newline.
	template<typename Type, typename ...Args>
	inline void println(const Type& curr, Args... args) {
		std::cout << curr << " ";
		println(args...);
	}

	/// Print the given argument to an output stream
	/// followed by a newline.
	template<typename Type>
	inline void fprintln(std::ostream& out, const Type& curr) {
		out << curr << "\n";
	}

	/// Print the given arguments to an output stream
	/// separated by a space and followed by a newline.
	template<typename Type, typename ...Args>
	inline void fprintln(std::ostream& out, const Type& curr, Args... args) {
		out << curr << " ";
		fprintln(out, args...);
	}


	/// Insert a data set of the given type from a stream,
	/// reading line by line until a line is equal to
	/// the terminator and parsing each line using the
	/// given function, returning the list of values.
	template<typename Type = real>
	inline vec<Type> readln(
		std::istream& in,
		const std::string& terminator,
		std::function<Type(std::string)> parse) {

		vec<Type> data;
		std::string line;
		Type value;

		while(true) {

			std::getline(in, line);

			// Stop reading when the terminator is reached.
			if(line == terminator)
				break;

			// Skip empty lines.
			if(line == "")
				continue;

			// Try to parse the line and add it to the vector.
			try {
				value = parse(line);
			} catch(...) {
				std::cout << "Input conversion error" << std::endl;
				continue;
			}

			data.push(value);
		}

		return data;
	}


	/// Insert a data set of the given type from a stream,
	/// reading line by line until a line is equal to
	/// the terminator and parsing each line as a real value,
	/// returning the list of values.
	inline vec<real> readln(
		std::istream& in, const std::string& terminator = "") {

		return readln<real>(in, terminator, [](std::string line) {
			std::replace(line.begin(), line.end(), ',', '.');
			return std::stod(line);
		});
	}


	/// Insert a data set of the given type from standard input,
	/// reading line by line until a line is equal to
	/// the terminator and parsing each line as a real value,
	/// returning the list of values.
	inline vec<real> readln(const std::string& terminator = "") {
		return readln(std::cin, terminator);
	}

}

#endif
