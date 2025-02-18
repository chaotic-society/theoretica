
# Build example programs

file(GLOB EXAMPLE_SOURCES *.cpp)

foreach(EXAMPLE_SRC IN LISTS EXAMPLE_SOURCES)

    get_filename_component(EXAMPLE_NAME ${EXAMPLE_SRC} NAME_WE)
    
    # Compile example program
    add_executable(${EXAMPLE_NAME} ${EXAMPLE_SRC})
    target_link_libraries(${EXAMPLE_NAME} theoretica)

endforeach()

# Add "examples" target to build all examples
add_custom_target(examples
    DEPENDS ${EXAMPLE_SOURCES}
    COMMENT "Building example programs..."
)
