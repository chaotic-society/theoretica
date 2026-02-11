
# Build test programs

file(GLOB TEST_SOURCES *.cpp)

# Include header files
target_include_directories(theoretica INTERFACE 
    ${CMAKE_CURRENT_SOURCE_DIR}/../src/
    ${CMAKE_CURRENT_SOURCE_DIR}/../chebyshev/src/
)

foreach(TEST_SRC IN LISTS TEST_SOURCES)

    get_filename_component(TEST_NAME ${TEST_SRC} NAME_WE)
    
    # Build test program
    add_executable(${TEST_NAME} ${TEST_SRC})
    target_link_libraries(${TEST_NAME} theoretica)
    
    # Add test executable to CTest
    add_test(
        NAME ${TEST_NAME}
        COMMAND ${TEST_NAME}
    )
    
    # Set passing criteria
    set_tests_properties(${TEST_NAME} PROPERTIES 
        WILL_FAIL OFF
        PASS_REGULAR_EXPRESSION "PASS|0 failed"
        FAIL_REGULAR_EXPRESSION "FAIL"
    )

endforeach()
