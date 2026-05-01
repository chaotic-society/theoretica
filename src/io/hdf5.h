
///
/// @file format_hdf5.h HDF5 file IO.
/// @note Because of its complexity and dependencies,
/// this module is not included automatically in the 'theoretica.h' header.
///

#ifndef THEORETICA_IO_FORMAT_HDF5_H
#define THEORETICA_IO_FORMAT_HDF5_H

#include <hdf5.h>

#include <string>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <cstring>
#include <sstream>
#include <algorithm>

#include "../core/constants.h"
#include "../algebra/vec.h"
#include "../algebra/mat.h"


namespace theoretica {
namespace io {


	// Data structures and representations

	/// Type of node inside an HDF5 file
	enum class HDF5NodeType {
		
		/// An unsupported node type
		UNKNOWN,

		/// A structural grouping containing other groups or datasets
		GROUP,

		/// A node containing multi-dimensional array data
		DATASET
	};


    /// @class hdf5_node
	/// Represents a single node (group or dataset) in the HDF5 file hierarchy.
	struct hdf5_node {

		/// The local name of the node (e.g. "my_dataset")
		std::string name;

		/// The absolute internal path to the node
		/// (e.g. "/group/subgroup/my_dataset")
		std::string path;

		/// The type of the node
		HDF5NodeType type;
		
		/// For datasets, the dimensions of the array
		std::vector<size_t> dimensions;

		/// List of metadata attribute names attached to the node
		std::vector<std::string> attributes;

		/// For groups, the child nodes indexed by their local name
		std::unordered_map<std::string, hdf5_node> children;


		/// Check whether the node is a group
		inline bool is_group() const noexcept {
			return type == HDF5NodeType::GROUP;
		}


		/// Check whether the node is a dataset
		inline bool is_dataset() const noexcept {
			return type == HDF5NodeType::DATASET;
		}


		/// Dictionary-like access to child nodes by name
		///
		/// @param child_name The local name of the child node
		/// @return Reference to the child node
		inline hdf5_node& operator[](const std::string& child_name) {
			return children[child_name];
		}


		/// Read-only dictionary-like access to child nodes by name
		///
		/// @param child_name The local name of the child node
		/// @return Const reference to the child node
		inline const hdf5_node& operator[](const std::string& child_name) const {
			return children.at(child_name);
		}
	};


	/// @class hdf5_handle
    /// RAII wrapper for managing HDF5 C-style handles, allowing direct API usage.
	struct hdf5_handle {

		/// Underlying ID obtained by opening
        /// a file or group.
		hid_t id;
		

		/// Initialize the handle with a given HDF5 ID,
        /// for example obtained by opening a file or group,
        /// defaulting to an invalid handle .
		hdf5_handle(hid_t id = H5I_INVALID_HID) : id(id) {}
		

		/// Safely decrements reference count and closes the handle if valid
		~hdf5_handle() {

			if (id >= 0 && H5Iis_valid(id))
				H5Idec_ref(id);
		}
		

		/// Prevent copying
		hdf5_handle(const hdf5_handle&) = delete;
		hdf5_handle& operator=(const hdf5_handle&) = delete;
		

		/// Move constructor
		hdf5_handle(hdf5_handle&& other) noexcept : id(other.id) {
			other.id = H5I_INVALID_HID;
		}


		/// Move assignment
		hdf5_handle& operator=(hdf5_handle&& other) noexcept {

			if (this != &other) {
				
				if (id >= 0 && H5Iis_valid(id))
					H5Idec_ref(id);
				
				id = other.id;
				other.id = H5I_INVALID_HID;
			}

			return *this;
		}


		/// Implicit conversion to hid_t
		operator hid_t() const {
			return id;
		}
	};


    namespace _internal {

        /// Helper template to map C++ types to HDF5 Native Types
        /// @tparam Type The C++ scalar type
        /// @return The equivalent HDF5 hid_t native type
        template<typename Type>
        inline hid_t hdf5_type() {
            throw std::invalid_argument("Unsupported type for HDF5 IO");
        }

        template<> inline hid_t hdf5_type<double>() { return H5T_NATIVE_DOUBLE; }
        template<> inline hid_t hdf5_type<float>()  { return H5T_NATIVE_FLOAT; }
        template<> inline hid_t hdf5_type<int>()    { return H5T_NATIVE_INT; }
        template<> inline hid_t hdf5_type<long>()   { return H5T_NATIVE_LONG; }
        template<> inline hid_t hdf5_type<unsigned int>() { return H5T_NATIVE_UINT; }


        /// Suppress HDF5 error messages by setting a custom dummy error handler.
        inline void suppress_errors() {
            H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
        }


		/// Removes a link (dataset or group) if it exists at the given path
		/// @param id The active file or group handle
		/// @param path The path to check and remove
		inline void remove_link(const hdf5_handle& id, const std::string& path) {

			suppress_errors();
			if (H5Lexists(id, path.c_str(), H5P_DEFAULT) > 0) {
				H5Ldelete(id, path.c_str(), H5P_DEFAULT);
			}
		}


		// Attributes

		// Count attributes callback for H5Aiterate2
		static herr_t attribute_callback(hid_t id, const char* attr_name, const H5A_info_t* info_ptr, void* op_data) {

			auto* attrs = (std::vector<std::string>*) op_data;
			attrs->emplace_back(attr_name);
			return 0;
		}


		// Load all attributes for a given object into a vector of strings
		inline void load_attributes(hid_t obj_id, std::vector<std::string>& attributes) {
			H5Aiterate2(obj_id, H5_INDEX_CRT_ORDER, H5_ITER_INC, NULL, attribute_callback, &attributes);
		}


		// Read and write attributes of various types, with error handling
		template<typename Type>
		inline void read_attribute(hid_t attr_id, Type& value) {

			if (H5Aread(attr_id, hdf5_type<Type>(), &value) < 0)
				throw std::runtime_error("Failed to read attribute");
		}


		// Read a string attribute
		inline void read_attribute(hid_t attr_id, std::string& value) {

			hdf5_handle type_id = H5Aget_type(attr_id);
			if (H5Tget_class(type_id) != H5T_STRING)
				throw std::runtime_error("Attribute is not a string");
			
			const int is_varstring = H5Tis_variable_str(type_id);
			if (is_varstring > 0) {

				char* buf = nullptr;
				hdf5_handle mem_type = H5Tcopy(H5T_C_S1);
				H5Tset_size(mem_type, H5T_VARIABLE);

				if (H5Aread(attr_id, mem_type, &buf) < 0)
					throw std::runtime_error("Failed to read var-string attribute");
				
				value = buf ? std::string(buf) : std::string();
				if (buf)
					H5free_memory(buf);

			} else if (is_varstring == 0) {

				const size_t size = H5Tget_size(type_id);
				std::vector<char> buf (size + 1, '\0');

				hdf5_handle mem_type = H5Tcopy(H5T_C_S1);
				H5Tset_size(mem_type, size);
				if (H5Aread(attr_id, mem_type, buf.data()) < 0)
					throw std::runtime_error("Failed to read string attribute");
				
				value = std::string(buf.data());

			} else {
				throw std::runtime_error("Failed to determine string attribute type");
			}
		}


		// Write an attribute of a generic type, creating it if it doesn't exist
		template<typename Type>
		inline void write_attribute(hid_t obj_id, const std::string& name, const Type& value) {

			hdf5_handle space_id = H5Screate(H5S_SCALAR);
			hdf5_handle attr_id = H5Acreate2(obj_id, name.c_str(), hdf5_type<Type>(), space_id, H5P_DEFAULT, H5P_DEFAULT);

			if (attr_id < 0)
				throw std::runtime_error("Failed to create attribute");

			if (H5Awrite(attr_id, hdf5_type<Type>(), &value) < 0)
				throw std::runtime_error("Failed to write attribute");
		}


		// Write a string attribute
		inline void write_attribute(hid_t obj_id, const std::string& name, const std::string& value) {

			hdf5_handle space_id = H5Screate(H5S_SCALAR);
			hdf5_handle type_id = H5Tcopy(H5T_C_S1);
			H5Tset_size(type_id, value.empty() ? 1 : value.length());
			H5Tset_strpad(type_id, H5T_STR_NULLTERM);
			
			hdf5_handle attr_id = H5Acreate2(obj_id, name.c_str(), type_id, space_id, H5P_DEFAULT, H5P_DEFAULT);
			if (attr_id < 0)
				throw std::runtime_error("Failed to create string attribute");

			if (H5Awrite(attr_id, type_id, value.c_str()) < 0)
				throw std::runtime_error("Failed to write string attribute");
		}


		// Write a C-style string attribute
		inline void write_attribute(hid_t obj_id, const std::string& name, const char* str) {
			write_attribute(obj_id, name, std::string(str));
		}

		
		// Recursively iterates over group members to build the hdf5_node tree structure
		inline herr_t iter_callback(hid_t id, const char* name, const H5L_info_t* info, void* opdata) {

			hdf5_node* parent = (hdf5_node*) opdata;
			hdf5_node child;
			child.name = name;

			// Build full child path
			if (parent->path == "/")
				child.path = "/" + std::string(name);
			else
				child.path = parent->path + "/" + std::string(name);

			// Open object by absolute child path from current location
			hdf5_handle obj_id = H5Oopen(id, name, H5P_DEFAULT);
			if (obj_id < 0) {
				child.type = HDF5NodeType::UNKNOWN;
				parent->children.emplace(child.name, std::move(child));
				return 0;
			}

			const H5I_type_t obj_type = H5Iget_type(obj_id);
			if (obj_type == H5I_GROUP) {

				child.type = HDF5NodeType::GROUP;

				// Read group attributes
				load_attributes(obj_id, child.attributes);

				// Recursion by opening the group and iterating inside it
				hdf5_handle gid = H5Gopen2(id, name, H5P_DEFAULT);
				if (gid >= 0)
					H5Literate(gid, H5_INDEX_NAME, H5_ITER_INC, NULL, iter_callback, &child);

			} else if (obj_type == H5I_DATASET) {

				child.type = HDF5NodeType::DATASET;

				// Dataset shape
				hdf5_handle space_id = H5Dget_space(obj_id);
				const int ndims = H5Sget_simple_extent_ndims(space_id);
				if (ndims > 0) {
					std::vector<hsize_t> dims (static_cast<size_t>(ndims));
					H5Sget_simple_extent_dims(space_id, dims.data(), NULL);
					child.dimensions = dims;
				}

				// Dataset attributes
				load_attributes(obj_id, child.attributes);

			} else {

				child.type = HDF5NodeType::UNKNOWN;
				load_attributes(obj_id, child.attributes);
			}

			parent->children.emplace(child.name, std::move(child));
			return 0;
		}


		// Helper function to recursively build a formatted
		// tree string representation of the HDF5 structure
		inline void build_tree(
			std::ostringstream& oss, const hdf5_node& node,
			const std::string& prefix, bool is_last, bool is_root) {

			// Current node line
			if (is_root)
				oss << (node.name.empty() ? "/" : node.name);
			else
				oss << prefix << (is_last ? "└── " : "├── ") << node.name;


			// Print dataset
			if (node.is_dataset()) {

				oss << " [";
				for (size_t i = 0; i < node.dimensions.size(); ++i) {
					oss << node.dimensions[i];
					if (i + 1 < node.dimensions.size())
						oss << ", ";
				}
				oss << "]";
			}


			// Show attribute count
			if (!node.attributes.empty()) {
				oss << " (@" << node.attributes.size() << " attrs)";
			}
			oss << "\n";


			if (node.children.empty())
				return;

			// Print group children recursively
			size_t i = 0;
			const std::string child_prefix = is_root ? "" : (prefix + (is_last ? "    " : "│   "));
			for (auto& pair : node.children) {

				const bool is_last_child = (node.children.size() == i + 1);
				build_tree(
					oss,
					pair.second,
					child_prefix,
					is_last_child,
					false
				);

				i++;
			}
		}
	}


	// Output formatting

	/// Generates a string representation of the HDF5 tree structure.
	///
	/// @param node The node to stringify (along its children)
	/// @return Formatted tree string
	inline std::string to_string(const hdf5_node& node) {
		std::ostringstream oss;
		_internal::build_tree(oss, node, "", true, true);
		return oss.str();
	}


	/// Prints the HDF5 file tree to a stream.
	///
	/// @param os The output stream
	/// @param node The target node
	/// @return Reference to the output stream
	inline std::ostream& operator<<(std::ostream& os, const hdf5_node& node) {
		return os << to_string(node);
	}


	// Public API using handles, for continued operations on an open file or group

	/// Open an HDF5 file with the given filename, returning a file handle.
	///
	/// @param filename The path to the HDF5 file
	/// @param write Whether to open the file with write permissions,
	/// creating it if it doesn't exist (defaults to false).
    inline hdf5_handle hdf5_open(const std::string& filename, bool write = false) {

        _internal::suppress_errors();
        hid_t file_id;

		// Open the file with write permissions,
		// creating it if it doesn't exist.
        if (write) {

            file_id = H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

            if (file_id < 0) {

                file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                if (file_id < 0)
                    throw std::runtime_error("Could not open or create HDF5 file: " + filename);
            }

        } else {

			// Open with read-only access.
            file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

            if (file_id < 0)
                throw std::runtime_error("Could not open HDF5 file: " + filename);
        }

        return hdf5_handle(file_id);
    }


	/// Recursively loads the structure of an active HDF5 group or file
	///
	/// @param id The open HDF5 handle representing the search root
	/// @return An hdf5_node tree hierarchy representing the file
	inline hdf5_node hdf5_load(const hdf5_handle& id) {

		_internal::suppress_errors();

		hdf5_node root;
		root.name = "/";
		root.path = "/";
		root.type = HDF5NodeType::GROUP;
		root.dimensions.clear();
		root.attributes.clear();
		root.children.clear();

		_internal::load_attributes(id, root.attributes);
		H5Literate(id, H5_INDEX_NAME, H5_ITER_INC, NULL, _internal::iter_callback, &root);

		return root;
	}


	/// Create a new group at the given path under an already open HDF5 location.
	///
	/// @param id Open HDF5 file/group identifier.
	/// @param path Absolute or relative group path to create.
	inline void hdf5_create_group(const hdf5_handle& id, const std::string& path) {

		_internal::suppress_errors();

		// Already exists, do nothing
		if (H5Lexists(id, path.c_str(), H5P_DEFAULT) > 0)
			return;

		hid_t gid = H5Gcreate2(id, path.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (gid < 0)
			throw std::runtime_error("Failed to create group: " + path);

		hdf5_handle h_group(gid);
	}


	/// Delete a group (link) at the given path under an already open HDF5 location.
	///
	/// @param id Open HDF5 file/group identifier.
	/// @param path Absolute or relative group path to delete.
	inline void hdf5_delete_group(const hdf5_handle& id, const std::string& path) {

		_internal::suppress_errors();

		// Do nothing if it doesn't exist
		if (H5Lexists(id, path.c_str(), H5P_DEFAULT) <= 0)
			return;

		if (H5Ldelete(id, path.c_str(), H5P_DEFAULT) < 0)
			throw std::runtime_error("Failed to delete group: " + path);
	}


	// Metadata handling

	/// Reads an attribute attached to a specific node.
	///
	/// @tparam Type The expected scalar type or std::string
	/// @param id The active file handle
	/// @param path Internal path to the node
	/// @param attr_name Name of the attribute to read
	/// @return The value of the attribute
	template<typename Type>
	inline Type hdf5_read_attribute(const hdf5_handle& id, const std::string& path, const std::string& attr_name) {

		_internal::suppress_errors();

		hid_t obj_id = H5Oopen(id, path.c_str(), H5P_DEFAULT);
		if (obj_id < 0)
			throw std::runtime_error("Cannot open node: " + path);

		hdf5_handle h_obj (obj_id);
		hid_t attr_id = H5Aopen(obj_id, attr_name.c_str(), H5P_DEFAULT);
		if (attr_id < 0)
			throw std::runtime_error("Cannot open attribute: " + attr_name);
		hdf5_handle h_attr (attr_id);
		
		Type value;
		_internal::read_attribute(attr_id, value);
		return value;
	}


	/// Deletes an attribute attached to a specific node, if it exists.
	///
	/// @param id The active file handle
	/// @param path Internal path to the node
	/// @param attr_name Name of the attribute to delete
	inline void hdf5_delete_attribute(const hdf5_handle& id, const std::string& path, const std::string& attr_name) {

		_internal::suppress_errors();

		hid_t obj_id = H5Oopen(id, path.c_str(), H5P_DEFAULT);
		if (obj_id < 0)
			throw std::runtime_error("Cannot open node: " + path);

		hdf5_handle h_obj (obj_id);

		if (H5Aexists(obj_id, attr_name.c_str()) > 0) {
			if (H5Adelete(obj_id, attr_name.c_str()) < 0)
				throw std::runtime_error("Failed to delete attribute: " + attr_name);
		}
	}


	/// Writes or overwrites an attribute attached to a specific node
	///
	/// @tparam Type The scalar type or std::string of the metadata
	/// @param id The active file handle
	/// @param path Internal path to the node
	/// @param attr_name Name of the attribute
	/// @param value The value to write
	template<typename Type>
	inline void hdf5_write_attribute(const hdf5_handle& id, const std::string& path, const std::string& attr_name, const Type& value) {

		_internal::suppress_errors();
		hid_t obj_id = H5Oopen(id, path.c_str(), H5P_DEFAULT);
		if (obj_id < 0)
			throw std::runtime_error("Cannot open node for metadata: " + path);

		hdf5_handle h_obj (obj_id);
		if (H5Aexists(obj_id, attr_name.c_str()) > 0)
			H5Adelete(obj_id, attr_name.c_str());

		_internal::write_attribute(obj_id, attr_name, value);
	}


	// Dataset IO for vectors and matrices

	/// Loads a 1D dataset array into a vector
	///
	/// @tparam Vector The vector container type (defaults to vec<real>).
	/// @param id The active file handle
	/// @param path Internal path to the dataset
	/// @return The populated vector
	template<typename Vector = vec<real>>
	inline Vector hdf5_read_vec(const hdf5_handle& id, const std::string& path) {

		using Type = vector_element_t<Vector>;
		
		_internal::suppress_errors();
		hid_t data_id = H5Dopen2(id, path.c_str(), H5P_DEFAULT);
		if (data_id < 0)
			throw std::runtime_error("Cannot open dataset: " + path);
		hdf5_handle h_dset (data_id);
		
		hdf5_handle h_space = H5Dget_space(data_id);
		if (H5Sget_simple_extent_ndims(h_space) != 1) 
			throw std::runtime_error("HDF5 Error: Node '" + path + "' is not a 1D dataset");
			
		hsize_t dims[1];
		H5Sget_simple_extent_dims(h_space, dims, NULL);
		
		Vector v (dims[0]);
		if (H5Dread(data_id, _internal::hdf5_type<Type>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
			throw std::runtime_error("Failed to read vector dataset");

		return v;
	}


	/// Writes a 1D vector to an HDF5 dataset, overwriting if it exists
	///
	/// @tparam Vector The vector container type.
	/// @param id The active file handle
	/// @param path Internal path to place the dataset
	/// @param v The vector data
	template<typename Vector>
	inline void hdf5_write_vec(const hdf5_handle& id, const std::string& path, const Vector& v) {

		using Type = vector_element_t<Vector>;
		_internal::remove_link(id, path);
		
		hsize_t dims[1] = { (hsize_t) v.size() };
		hdf5_handle h_space = H5Screate_simple(1, dims, NULL);
		
		hid_t data_id = H5Dcreate2(id, path.c_str(), _internal::hdf5_type<Type>(), h_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (data_id < 0)
			throw std::runtime_error("Failed to create dataset: " + path);
		hdf5_handle h_dset(data_id);
		
		if (H5Dwrite(data_id, _internal::hdf5_type<Type>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, v.data()) < 0)
			throw std::runtime_error("Failed to write vector dataset");
	}


	/// Deletes a dataset at the given path if it exists.
	///
	/// @param id The active file handle
	/// @param path Internal path to the dataset
	inline void hdf5_delete_dataset(const hdf5_handle& id, const std::string& path) {

		_internal::suppress_errors();

		if (H5Lexists(id, path.c_str(), H5P_DEFAULT) > 0) {
			if (H5Ldelete(id, path.c_str(), H5P_DEFAULT) < 0)
				throw std::runtime_error("Failed to delete dataset: " + path);
		}
	}


	/// Loads a 2D dataset array into a matrix
	///
	/// @tparam Matrix The matrix container type (defaults to mat<real>).
	/// @param id The active file handle
	/// @param path Internal path to the dataset
	/// @return The populated matrix
	template<typename Matrix = mat<real>>
	inline Matrix hdf5_read_mat(const hdf5_handle& id, const std::string& path) {

		using Type = matrix_element_t<Matrix>;
		_internal::suppress_errors();

		hid_t data_id = H5Dopen2(id, path.c_str(), H5P_DEFAULT);
		if (data_id < 0)
			throw std::runtime_error("Cannot open dataset: " + path);
		hdf5_handle h_dset(data_id);
		
		hdf5_handle h_space = H5Dget_space(data_id);
		if (H5Sget_simple_extent_ndims(h_space) != 2) 
			throw std::runtime_error("HDF5 Error: Node '" + path + "' is not a 2D dataset");
			
		hsize_t dims[2];
		H5Sget_simple_extent_dims(h_space, dims, NULL);
		
		Matrix m (dims[0], dims[1]);
		if (H5Dread(data_id, _internal::hdf5_type<Type>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, m.data()) < 0)
			throw std::runtime_error("Failed to read matrix dataset");

		return m;
	}


	/// Writes a 2D matrix to an HDF5 dataset, overwriting if it exists
	///
	/// @tparam Matrix A generic matrix type. Must support .rows(), .cols(), and .data().
	/// @param id The active file handle
	/// @param path Internal path to place the dataset
	/// @param m The matrix data
	template<typename Matrix>
	inline void hdf5_write_mat(const hdf5_handle& id, const std::string& path, const Matrix& m) {

		using Type = matrix_element_t<Matrix>;
		_internal::remove_link(id, path);
		
		hsize_t dims[2] = { hsize_t(m.rows()), hsize_t(m.cols()) };
		hdf5_handle h_space = H5Screate_simple(2, dims, NULL);
		
		hid_t data_id = H5Dcreate2(id, path.c_str(), _internal::hdf5_type<Type>(), h_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
		if (data_id < 0)
			throw std::runtime_error("Failed to create dataset: " + path);
		
		hdf5_handle h_dset (data_id);
		
		if (H5Dwrite(data_id, _internal::hdf5_type<Type>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, m.data()) < 0) {
			throw std::runtime_error("Failed to write matrix dataset");
		}
	}


	/// @class hdf5_file
	/// High-level interface for managing an HDF5 file, its structure, and operations.
	class hdf5_file {
	private:
		
		/// Filename for reference and error messages
		std::string m_filename;

		/// Internal handle for the open HDF5 file
		hdf5_handle m_file_id; 

		/// Root node representing the entire file structure
		hdf5_node m_root;

	public:

		/// Open the file persistently
		///
		/// @param filename Path to the target HDF5 file
		/// @param write If true, opens the file with write permissions (creating it if it doesn't exist)
		explicit hdf5_file(const std::string& filename, bool write = false) : m_filename(filename) {
			m_file_id = hdf5_open(m_filename, write);
			refresh();
		}


		/// Refresh the cached tree structure from disk
		void refresh() {

			// Always rebuild from scratch
			m_root = hdf5_node{};
			m_root.name = "/";
			m_root.path = "/";
			m_root.type = HDF5NodeType::GROUP;
			m_root.dimensions.clear();
			m_root.attributes.clear();
			m_root.children.clear();

			// Rebuild root from open file handle
			m_root = hdf5_load(m_file_id);

			// Keep root name as filename
			m_root.name = m_filename;
		}


		/// Close the HDF5 file immediately
		void close() {

			H5Fclose(m_file_id);
			m_file_id = H5I_INVALID_HID;

			// Clear cached structure
			m_root = hdf5_node();
		}


		/// Get the name of the file
		/// @return String reference of the path
		const std::string& filename() const {
			return m_filename;
		}


		/// Get the base node containing the entire loaded file structure
		/// @return A const hdf5_node reference representing the root "/"
		const hdf5_node& root() const {
			return m_root;
		}


		/// Get the raw ID of the HDF5 file handle for direct API access
		/// @return Internal HDF5 file identifier
		hid_t id() const {
			return m_file_id;
		}


		/// Access children of the root node by reference using dictionary-like syntax.
		/// @param child_name Node name
		hdf5_node& operator[](const std::string& child_name) {
			return m_root[child_name];
		}


		/// Access children of the root node by const reference using dictionary-like syntax.
		/// @param child_name Node name
		const hdf5_node& operator[](const std::string& child_name) const {
			return m_root[child_name];
		}


		/// Create a new group inside the file.
		/// @param path Absolute or relative group path.
		inline void create_group(const std::string& path) {
			hdf5_create_group(m_file_id, path);
		}


		/// Create a new group using a parent node and a child name.
		/// @param parent Parent node under which the new group will be created.
		/// @param child_name Name of the child group.
		inline void create_group(const hdf5_node& parent, const std::string& child_name) {

			const std::string path = (parent.path == "/")
				? ("/" + child_name)
				: (parent.path + "/" + child_name);

			hdf5_create_group(m_file_id, path);
		}


		/// Delete a group inside the file.
		/// @param path Absolute or relative group path.
		inline void delete_group(const std::string& path) {
			hdf5_delete_group(m_file_id, path);
		}


		/// Delete a group referenced by a node.
		/// @param node Node to delete.
		inline void delete_group(const hdf5_node& node) {
			hdf5_delete_group(m_file_id, node.path);
		}


		// Vector operations

		/// Loads a 1D dataset array into a vector from the given path
		/// inside the HDF5 file.
		///
		/// @tparam Vector The return container type (defaults to vec<real>).
		/// @param path Internal path to the dataset
		/// @return The vector populated with the dataset values
		template<typename Vector = vec<real>>
		Vector read_vec(const std::string& path) const {
			return hdf5_read_vec<Vector>(m_file_id, path);
		}


		/// Loads a 1D dataset array into a vector from the given node
		/// inside the HDF5 file.
		///
		/// @tparam Vector The return container type (defaults to vec<real>).
		/// @param node Node representing the dataset to read
		/// @return The vector populated with the dataset values
		template<typename Vector = vec<real>>
		Vector read_vec(const hdf5_node& node) const {
			return hdf5_read_vec<Vector>(m_file_id, node.path);
		}


		/// Writes a 1D vector to an HDF5 dataset at the given path, overwriting if it exists.
		///
		/// @tparam Vector The generic container type.
		/// @param path Internal path to place the dataset
		/// @param v The vector data to write
		template<typename Vector>
		void write_vec(const std::string& path, const Vector& v) {
			hdf5_write_vec(m_file_id, path, v);
		}

		/// Writes a 1D vector to the HDF5 dataset represented by the given node, overwriting if it exists.
		///
		/// @tparam Vector The generic container type.
		/// @param node HDF5 node identifying the dataset to write
		/// @param v The vector data to write
		template<typename Vector>
		void write_vec(const hdf5_node& node, const Vector& v) {
			hdf5_write_vec(m_file_id, node.path, v);
		}


		// Matrix operations
		
		/// Read a 2D dataset array into a matrix from the given path.
		///
		/// @tparam Matrix The return container type (defaults to mat<real>).
		/// @param path Internal path to the dataset
		/// @return The matrix populated with the dataset values
		template<typename Matrix = mat<real>>
		Matrix read_mat(const std::string& path) const {
			return hdf5_read_mat<Matrix>(m_file_id, path);
		}
		

		/// Read a 2D dataset array into a matrix from the given node.
		///
		/// @tparam Matrix The return container type (defaults to mat<real>).
		/// @param node Node representing the dataset to read
		/// @return The matrix populated with the dataset values
		template<typename Matrix = mat<real>>
		Matrix read_mat(const hdf5_node& node) const {
			return hdf5_read_mat<Matrix>(m_file_id, node.path);
		}
		

		/// Write a 2D matrix to an HDF5 dataset at the given path, overwriting if it exists.
		///
		/// @tparam Matrix The generic container type (deduced). Must support .rows(), .cols(), and .data().
		/// @param path Internal path to place the dataset
		/// @param m The matrix data to write
		template<typename Matrix>
		void write_mat(const std::string& path, const Matrix& m) {
			hdf5_write_mat(m_file_id, path, m);
		}
		

		/// Write a 2D matrix to an HDF5 dataset at the given node, overwriting if it exists.
		///
		/// @tparam Matrix The generic container type (deduced). Must support .rows(), .cols(), and .data().
		/// @param node Node representing the dataset to write
		/// @param m The matrix data to write
		template<typename Matrix>
		void write_mat(const hdf5_node& node, const Matrix& m) {
			hdf5_write_mat(m_file_id, node.path, m);
		}


		/// Delete a dataset at the given path if it exists.
		///
		/// @param path Internal path to the dataset
		void delete_dataset(const std::string& path) {
			hdf5_delete_dataset(m_file_id, path);
		}


		/// Delete a dataset at the given node if it exists.
		///
		/// @param node Node representing the dataset to delete
		void delete_dataset(const hdf5_node& node) {
			hdf5_delete_dataset(m_file_id, node.path);
		}


		// Metadata operations
		
		/// Read metadata (attribute) attached to a specific node by path
		///
		/// @tparam Type The expected scalar type or std::string of the attribute
		/// @param path Internal path to the node
		/// @param attr_name Name of the attribute to read
		/// @return The value of the attribute
		template<typename Type>
		Type read_attribute(const std::string& path, const std::string& attr_name) const {
			return hdf5_read_attribute<Type>(m_file_id, path, attr_name);
		}
		

		/// Read metadata (attribute) attached to a specific node by reference
		///
		/// @tparam Type The expected scalar type or std::string of the attribute
		/// @param node Node representing the object with the attribute
		/// @param attr_name Name of the attribute to read
		/// @return The value of the attribute
		template<typename Type>
		Type read_attribute(const hdf5_node& node, const std::string& attr_name) const {
			return hdf5_read_attribute<Type>(m_file_id, node.path, attr_name);
		}
		

		/// Write or overwrite metadata (attribute) attached to a specific node by path
		///
		/// @tparam Type The scalar type or std::string of the attribute
		/// @param path Internal path to the node
		/// @param attr_name Name of the attribute
		/// @param value The value to write
		template<typename Type>
		void write_attribute(const std::string& path, const std::string& attr_name, const Type& value) {
			hdf5_write_attribute<Type>(m_file_id, path, attr_name, value);
		}
		

		/// Write or overwrite metadata (attribute) attached to a specific node by reference
		///
		/// @tparam Type The scalar type or std::string of the attribute
		/// @param node Node representing the object with the attribute
		/// @param attr_name Name of the attribute
		/// @param value The value to write
		template<typename Type>
		void write_attribute(const hdf5_node& node, const std::string& attr_name, const Type& value) {
			hdf5_write_attribute<Type>(m_file_id, node.path, attr_name, value);
		}


		/// Delete an attribute attached to a specific node by path, if it exists.
		///
		/// @param path Internal path to the node
		/// @param attr_name Name of the attribute to delete
		void delete_attribute(const std::string& path, const std::string& attr_name) {
			hdf5_delete_attribute(m_file_id, path, attr_name);
		}


		/// Delete an attribute attached to a specific node by path, if it exists.
		///
		/// @param node Node representing the object with the attribute
		/// @param attr_name Name of the attribute to delete
		void delete_attribute(const hdf5_node& node, const std::string& attr_name) {
			hdf5_delete_attribute(m_file_id, node.path, attr_name);
		}


		/// Converts an entire HDF5 file structure into a formatted tree string.
		///
		/// @param file The opened HDF5 file to format
		/// @return Formatted string hierarchy
		inline std::string to_string() const {
			return io::to_string(m_root);
		}


		/// Print the HDF5 file structure to a stream in a formatted tree representation.
		///
		/// @param os Standard output stream
		/// @param file The HDF5 wrapper instance
		/// @return Reference to the stream
		inline friend std::ostream& operator<<(std::ostream& os, const hdf5_file& file) {
			return os << file.to_string();
		}
	};

}}

#endif
