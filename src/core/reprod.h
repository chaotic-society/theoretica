
///
/// @file reprod.h Reproducibility features for reliable scientific computing.
///

#ifndef THEORETICA_REPROD_H
#define THEORETICA_REPROD_H

#include <string>


namespace theoretica {

/// @namespace Reproducibility features.
namespace reprod {

    
    /// @class environment
    /// Structure containing information about the build environment,
    /// such as operating system and compiler, for reproducibility purposes.
    struct environment {
        
        // Operating system
        std::string os {""};

        // Architecture
        std::string arch {""};

        // Compiler
        std::string compiler {""};

        // Compiler version
        std::string compiler_version {""};

        // C++ Standard
        std::string cpp_standard {""};

        // Build date
        std::string build_date {""};

        // AVX2 support
        bool has_avx2 {false};

        // AVX512 support
        bool has_avx512 {false};

        // CUDA support
        bool has_cuda {false};

        // OpenMP support
        bool has_omp {false};


        /// Convert the environment to human-readable string representation.
        inline std::string to_string() const {

            std::string str = "Environment Information\n";
            str += "OS: " + os + "\n";
            str += "Architecture: " + arch + "\n";
            str += "Compiler: " + compiler + " " + compiler_version + "\n";
            str += "C++ Standard: " + cpp_standard + "\n";
            str += "Features: ";

            if (has_avx2) str += "AVX2 ";
            if (has_avx512) str += "AVX512 ";
            if (has_cuda) str += "CUDA ";
            if (has_omp) str += "OpenMP ";

            str += "\n";
            //str += "Theoretica Version:  " + lib_version + "\n";
            str += "Build Date: " + build_date + "\n";
            
            return str;
        }

        /// Convert the environment to string representation.
        inline operator std::string() {
            return to_string();
        }


        /// Stream the environment in string representation to an output stream (std::ostream)
        inline friend std::ostream& operator<<(std::ostream& out, const environment& obj) {
            return out << obj.to_string();
        }

    };
    

    /// Get an environment structure holding information
    /// about the current executing environment.
    inline environment get_env() {

        environment env;

        // Operating system
        #ifdef _WIN32
            env.os = "Windows";
            #ifdef _WIN64
                env. arch = "x86_64";
            #else
                env.arch = "x86";
            #endif
        #elif defined(__APPLE__)
            env.os = "macOS";
            #ifdef __arm64__
                env.arch = "arm64";
            #else
                env.arch = "x86_64";
            #endif
        #elif defined(__linux__)
            env.os = "Linux";
            #ifdef __x86_64__
                env.arch = "x86_64";
            #elif defined(__aarch64__)
                env.arch = "arm64";
            #else
                env.arch = "";
            #endif
        #endif
        
        // Compiler detection
        #ifdef __GNUC__
            env.compiler = "gcc";
            env.compiler_version = std::to_string(__GNUC__) + "."
                + std::to_string(__GNUC_MINOR__) + "."
                + std::to_string(__GNUC_PATCHLEVEL__);
        #elif defined(__clang__)
            env.compiler = "clang";
            env.compiler_version = std::to_string(__clang_major__) + "."
            + std:: to_string(__clang_minor__) + "."
            + std::to_string(__clang_patchlevel__);
        #elif defined(_MSC_VER)
            env.compiler = "msvc";
            env.compiler_version = std::to_string(_MSC_VER);
        #endif
        
        // C++ Standard
        #if __cplusplus == 202002L
            env.cpp_standard = "C++20";
        #elif __cplusplus == 201703L
            env.cpp_standard = "C++17";
        #elif __cplusplus == 201402L
            env.cpp_standard = "C++14";
        #else
            env.cpp_standard = "C++11 or earlier";
        #endif
        
        // CPU features
        #ifdef __AVX2__
            env.has_avx2 = true;
        #endif

        #ifdef __AVX512F__
            env.has_avx512 = true;
        #endif

        #ifdef _OPENMP
            env.has_omp = true;
        #endif
        
        // GPU features
        #ifdef __CUDA_ARCH__
            env.has_cuda = true;
        #endif
        
        // Theoretica versioning
        //env.lib_version = std::string(THEORETICA_MAJOR) + "." + std::string(THEORETICA_MAJOR) + "." + std::string(THEORETICA_MAJOR);
        env.build_date = __DATE__ " " __TIME__;
        
        return env;

    }

}}

#endif
