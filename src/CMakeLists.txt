file(GLOB FLOXER_SOURCE_FILES RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/*.cpp)
list(REMOVE_ITEM FLOXER_SOURCE_FILES main.cpp)

add_library ("${PROJECT_NAME}_lib" ${FLOXER_SOURCE_FILES})
target_link_libraries (
    "${PROJECT_NAME}_lib" PUBLIC
    # add dependecies here
    sharg::sharg
    ivio::ivio
    ivsigma::ivsigma
    fmindex-collection::fmindex-collection
    cereal::cereal
    spdlog::spdlog
)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC "../include")
target_compile_options ("${PROJECT_NAME}_lib" PUBLIC "-pedantic" "-Wall" "-Wextra")
if (NOT ${PROJECT_NAME}_NO_OPENMP)
    target_compile_options("${PROJECT_NAME}_lib" PUBLIC "-fopenmp")
endif ()

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries("${PROJECT_NAME}" PUBLIC "${PROJECT_NAME}_lib")
if (NOT ${PROJECT_NAME}_NO_OPENMP)
    target_link_libraries("${PROJECT_NAME}" PRIVATE "-fopenmp")
endif ()
