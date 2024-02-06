# An interface to reuse common properties.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE sharg::sharg)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE ivio::ivio)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE ivsigma::ivsigma)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE fmindex-collection::fmindex-collection)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "../include")
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

add_library (cli_lib STATIC cli.cpp)
target_link_libraries (cli_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (io_lib STATIC io.cpp)
target_link_libraries (io_lib PUBLIC "${PROJECT_NAME}_interface")

add_library (pex_lib STATIC pex.cpp)
target_link_libraries (pex_lib PUBLIC "${PROJECT_NAME}_interface")

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries (
    "${PROJECT_NAME}"
    PRIVATE
    cli_lib
    io_lib
    pex_lib
)
