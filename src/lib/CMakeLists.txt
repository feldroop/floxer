file(GLOB FLOXER_LIB_SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

add_library ("${PROJECT_NAME}_lib" ${FLOXER_LIB_SOURCE_FILES})

target_link_libraries (
    "${PROJECT_NAME}_lib" PUBLIC
    sharg::sharg
    spdlog::spdlog
    ivio::ivio
    ivsigma::ivsigma
    fmindex-collection::fmindex-collection
    cereal::cereal
    seqan3::seqan3
    interval-tree
    BS_thread_pool
)

target_include_directories ("${PROJECT_NAME}_lib" PUBLIC "../../include")

target_compile_options ("${PROJECT_NAME}_lib" PUBLIC "-pedantic" "-Wall" "-Wextra")
