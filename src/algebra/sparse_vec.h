
///
/// @file sparse_vec.h Sparse vector class implemented as a hashmap
///

#ifndef THEORETICA_SPARSE_VEC_H
#define THEORETICA_SPARSE_VEC_H

#include <map>
#include "../core/constants.h"


namespace theoretica {

	/// @class sparse_vec Sparse vector class implemented as a hashmap
	template<typename T = real, typename IndexType = size_t>
	class sparse_vec {
		public:

			/// Data is stored in a hashmap as (index, value)
			std::map<IndexType, T> data;


			/// Initialize as the zero vector
			sparse_vec() {}
			

			/// Initialize the sparse vector from a list
			/// of pairs (index, value)
			template<typename Pair>
			sparse_vec(const std::initializer_list<Pair>& pairs) {

				for (auto& pair : pairs)
					data[pair[0]] = pair[1];
			}


			~sparse_vec() = default;


			/// Add two sparse vectors
			inline sparse_vec<T> operator+(const sparse_vec<T>& other) {

				sparse_vec res;
				res.data = this->data;

				for (auto& pair : other)
					res[pair.first] += pair.second;

				return res;
			}


			/// Dot product of two sparse vectors
			inline T operator*(const sparse_vec<T>& other) {

				T res = T(0.0);

				for (auto& pair : data)
					res += pair.second * other[pair.first];

				return res;
			}


			/// Get the n-th element
			inline T operator[](IndexType index) const {

				const auto iter = data.find(index);

				if(iter == data.end())
					return 0;

				return iter->second;
			}


			/// Access the n-th element
			inline T& operator[](IndexType index) {
				return data[index];
			}


#ifndef THEORETICA_NO_PRINT

		/// Convert the vector to string representation
		inline std::string to_string(
			const std::string& separator = ", ",
			bool parenthesis = true) const {

			std::stringstream res;

			// if(parenthesis)
			// 	res << "(";

			for (auto& pair : this->data) {
				// res << pair.second;
				// res << separator;

				res << "(" << pair.first << ", "
					<< pair.second << ") ";

			}

			// if(parenthesis)
			// 	res << ")";

			return res.str();
		}


		/// Stream the vector in string representation to an output stream (std::ostream)
		inline friend std::ostream& operator<<(
			std::ostream& out, const sparse_vec<T>& obj) {
			return out << obj.to_string();
		}

#endif

	};

}

#endif
