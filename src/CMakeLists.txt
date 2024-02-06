# An interface to reuse common properties.
add_library ("${PROJECT_NAME}_interface" INTERFACE)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE sharg::sharg)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE ivio::ivio)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE ivsigma::ivsigma)
target_link_libraries ("${PROJECT_NAME}_interface" INTERFACE fmindex-collection::fmindex-collection)
# target_include_directories ("${PROJECT_NAME}_interface" INTERFACE "${CMAKE_CURRENT_LIST_DIR}/../include")
target_compile_options ("${PROJECT_NAME}_interface" INTERFACE "-pedantic" "-Wall" "-Wextra")

#An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC floxer.cpp)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC "${PROJECT_NAME}_interface")

add_executable ("${PROJECT_NAME}" main.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
