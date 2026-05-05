
///
/// @file io/error.h Error handling for IO operations.
///

#ifndef THEORETICA_IO_ERROR_H
#define THEORETICA_IO_ERROR_H

#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)
#include <exception>
#ifndef THEORETICA_NO_PRINT
#include <sstream>
#endif
#endif

#include <string>
#include <cerrno>


namespace theoretica {
namespace io {


    /// IO error enumeration
    enum class IoError : int {

        /// No error
        None = 0x00,

        /// File or directory not found
        FileNotFound = 0x01,

        /// Permission denied when accessing the file or directory
        PermissionDenied = 0x02,

        /// Error occurred while reading from the file or stream
        ReadError = 0x04,

        /// Error occurred while writing to the file or stream
        WriteError = 0x08,

        /// The file format is invalid or the data is corrupted
        FormatError = 0x10,
        
        /// The end of the file was reached unexpectedly during a read operation
        EndOfFile = 0x20
    };


    /// Convert an IoError class enum to conventional errno codes.
    inline int to_errno(IoError err) {
        switch(err) {
            case IoError::None: return 0;
            case IoError::FileNotFound: return ENOENT;
            case IoError::PermissionDenied: return EACCES;
            case IoError::ReadError: return EIO;
            case IoError::WriteError: return EIO;
            case IoError::FormatError: return EINVAL;
            case IoError::EndOfFile: return ENODATA;
            default: return 0;
        }
    }


    /// Convert an IoError class enum to a string description.
    inline const char* to_cstring(IoError err) {
        switch (err) {
            case IoError::None: return "No error";
            case IoError::FileNotFound: return "File not found";
            case IoError::PermissionDenied: return "Permission denied";
            case IoError::ReadError: return "Failed to read data from IO stream";
            case IoError::WriteError: return "Failed to write data to IO stream";
            case IoError::FormatError: return "Invalid or corrupted file format";
            case IoError::EndOfFile: return "Unexpected end of file";
            default: return "Unknown IO error";
        }
    }

    inline std::string to_string(IoError err) {
        return to_cstring(err);
    }


#if defined(THEORETICA_THROW_EXCEPTIONS) || defined(THEORETICA_ONLY_EXCEPTIONS)

    class io_exception : public std::exception {
    private:

        IoError err;
        std::string func_name;
        std::string code_file_name;
        unsigned int code_line;
        std::string resource;

    public:

        io_exception(IoError a_err, const std::string& a_func_name,
            const std::string& a_code_file_name, unsigned int a_code_line, const std::string& a_target_file)
                : err(a_err), func_name(a_func_name), code_file_name(a_code_file_name),
                  code_line(a_code_line), resource(a_target_file) {}

        ~io_exception() = default;

        inline const char* what() const noexcept override {

            std::string str = code_file_name
                + "(" + std::to_string(code_line) + "):"
                + func_name + "(\"" + resource + "\"): "
                + to_cstring(err);

            return str.c_str();
        }


        /// Get the error code associated with the exception
        inline IoError err_code() const {
            return err;
        }


        /// Get the name of the throwing function
        inline std::string get_function_name() const {
            return func_name;
        }


        /// Get the name of the file in which the exception was thrown
        inline std::string get_file_name() const {
            return code_file_name;
        }


        /// Get the line number at which the exception was thrown
        inline unsigned int get_line_number() const {
            return code_line;
        }


        /// Get the resource (e.g., filename) associated with the exception
        inline std::string get_resource_id() const {
            return resource;
        }


#ifndef THEORETICA_NO_PRINT

        inline std::string to_string() const {
            std::stringstream err_str;
            err_str << code_file_name << "(" << code_line << "):";
            err_str << func_name << "(\"" << resource << "\"): " << to_cstring(err);
            return err_str.str();
        }

        inline operator std::string() const {
            return to_string();
        }

        inline friend std::ostream& operator<<(std::ostream& out, const io_exception& obj) {
            return out << obj.to_string();
        }
#endif

    };

#endif
}}


// Only throw exceptions, without modifying errno
#ifdef THEORETICA_ONLY_EXCEPTIONS

#define TH_IO_ERROR(F_NAME, RESOURCE_ID, EXCEPTION) \
    { throw io::io_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, RESOURCE_ID); }

// Throw exceptions and modify errno
#elif defined(THEORETICA_THROW_EXCEPTIONS)

#define TH_IO_ERROR(F_NAME, RESOURCE_ID, EXCEPTION) \
    { errno = io::to_errno(EXCEPTION); \
    throw io::io_exception(EXCEPTION, F_NAME, __FILE__, __LINE__, RESOURCE_ID); }

// Modify errno only by default
#else

#define TH_IO_ERROR(F_NAME, RESOURCE_ID, EXCEPTION) \
    { errno = io::to_errno(EXCEPTION); }

#endif


#endif