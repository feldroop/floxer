# this interface is used to give the same properties to all libraries
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries (
    "${PROJECT_NAME}_interface" INTERFACE
    # add dependecies here
    sharg::sharg
    ivio::ivio
    ivsigma::ivsigma
    fmindex-collection::fmindex-collection
    cereal::cereal
    spdlog::spdlog
)
target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "../include")
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")
if (NOT ${PROJECT_NAME}_NO_OPENMP)
    target_compile_options("${PROJECT_NAME}_interface" INTERFACE "-fopenmp")
endif ()

file(GLOB FLOXER_SOURCE_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable ("${PROJECT_NAME}" main.cpp)
if (NOT ${PROJECT_NAME}_NO_OPENMP)
    target_link_libraries("${PROJECT_NAME}" PRIVATE "-fopenmp")
endif ()

# every file in this folder, except for main.cpp is a static library that will
# be linked with the dependencies and main will be linked with it
foreach (_src_file ${FLOXER_SOURCE_FILES})
    get_filename_component (_lib_name ${_src_file} NAME_WE)
    if (NOT ("${_lib_name}" STREQUAL "main"))
        add_library ("${_lib_name}" STATIC "${_src_file}")
        target_link_libraries ("${_lib_name}" PUBLIC "${PROJECT_NAME}_interface")
        list(APPEND FLOXER_LIBRARIES "${_lib_name}")
    endif ()
endforeach ()

# this is actually necessary twice to make sure that all the links are found
# the linker lld works without this, but ld needs it
target_link_libraries("${PROJECT_NAME}" PRIVATE ${FLOXER_LIBRARIES})
target_link_libraries("${PROJECT_NAME}" PRIVATE ${FLOXER_LIBRARIES})
