# An interface to reuse common properties.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries (
    "${PROJECT_NAME}_interface" INTERFACE
    sharg::sharg
    ivio::ivio
    ivsigma::ivsigma
    fmindex-collection::fmindex-collection
    cereal::cereal
    spdlog::spdlog
)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE )
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE )
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE )
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE )
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE )
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "../include")
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra" "-fopenmp")

add_library (cli_lib STATIC cli.cpp)
target_link_libraries (cli_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (input_lib STATIC input.cpp)
target_link_libraries (input_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (output_lib STATIC output.cpp)
target_link_libraries (output_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (pex_lib STATIC pex.cpp)
target_link_libraries (pex_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (search_lib STATIC search.cpp)
target_link_libraries (search_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (alignment_lib STATIC alignment.cpp)
target_link_libraries (alignment_lib PUBLIC "${PROJECT_NAME}_interface")

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries (
    "${PROJECT_NAME}"
    PRIVATE
    cli_lib
    input_lib
    output_lib
    pex_lib
    search_lib
    alignment_lib
    -fopenmp
)
