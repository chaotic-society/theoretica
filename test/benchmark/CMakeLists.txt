
# Build benchmark programs

file(GLOB BENCHMARK_SOURCES *.cpp)

# Include header files
target_include_directories(theoretica INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/
    ${CMAKE_CURRENT_SOURCE_DIR}/../chebyshev/src/
)

foreach(BENCHMARK_SRC IN LISTS BENCHMARK_SOURCES)

    get_filename_component(BENCHMARK_NAME ${BENCHMARK_SRC} NAME_WE)
    
    # Build benchmark program
    add_executable(${BENCHMARK_NAME} ${BENCHMARK_SRC})
    target_link_libraries(${BENCHMARK_NAME} theoretica)

endforeach()

# Add "benchmark" target to run all benchmarks
add_custom_target(benchmark
    DEPENDS ${BENCHMARK_SOURCES}
    COMMENT "Building benchmarks..."
)
