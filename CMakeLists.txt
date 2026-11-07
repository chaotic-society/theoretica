
cmake_minimum_required(VERSION 3.12)
project(Theoretica)

# C++ Standard
set(CXX_STANDARD 20 CACHE STRING "C++ Standard (14,17,20)")
set(CMAKE_CXX_STANDARD ${CXX_STANDARD})
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add header-only library target
add_library(theoretica INTERFACE)
target_include_directories(theoretica 
    INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/src/
)

# Add examples
add_subdirectory(examples)

# Add tests
enable_testing()
add_subdirectory(test/prec)

# Add benchmarks
add_subdirectory(test/benchmark)
