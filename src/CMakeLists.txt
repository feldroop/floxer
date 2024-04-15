add_library ("${PROJECT_NAME}_interface" INTERFACE)
# add dependecies here
target_link_libraries (
    "${PROJECT_NAME}_interface" INTERFACE
    sharg::sharg
    ivio::ivio
    ivsigma::ivsigma
    fmindex-collection::fmindex-collection
    cereal::cereal
    spdlog::spdlog
)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "../include")
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra" "-fopenmp")

file(GLOB FLOXER_SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries("${PROJECT_NAME}" PRIVATE "-fopenmp")

# every file in this folder, except for main.cpp is a static library that will
# be linked to the dependencies and main will be linked to it
foreach (_src_file ${FLOXER_SOURCE_FILES})
    get_filename_component (_lib_name ${_src_file} NAME_WE)
    if (NOT ("${_lib_name}" STREQUAL "main"))
        add_library ("${_lib_name}" STATIC "${_src_file}")
        target_link_libraries ("${_lib_name}" PUBLIC "${PROJECT_NAME}_interface")
        target_link_libraries ("${PROJECT_NAME}" PRIVATE "${_lib_name}")
    endif ()
endforeach ()
