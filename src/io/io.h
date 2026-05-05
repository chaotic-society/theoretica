
///
/// @file io.h Functions for standard input and output.
///

#ifndef THEORETICA_IO_H
#define THEORETICA_IO_H

#include <iostream>
#include <string>


namespace theoretica {

	/// @namespace theoretica::io Input and output module.
	namespace io {

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


		/// Read a line from standard input, up to a line return.
		///
		/// @return The line as a string.
		inline std::string readln() {

			std::string str;
			std::getline(std::cin, str);

			return str;
		}


		/// Read a single object from standard input, ended by a line return.
		/// The object needs to support streaming from an std::istream.
		template<typename Type>
		inline void readln(Type& last) {
			std::cin >> last;
		}

		/// Read objects from standard input, ended by a line return.
		/// The object needs to support streaming from an std::istream.
		///
		/// For example, calling readln(x1, x2, x3) reads a line of the form:
		/// "1.0 2.0 3.0\n" into variables x1, x2 and x3.
		template<typename Type, typename ...Args>
		inline void readln(Type& curr, Args&... args) {
			std::cin >> curr;
			readln(args...);
		}

	}
}


#endif
