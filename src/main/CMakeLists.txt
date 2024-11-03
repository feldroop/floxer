file(GLOB FLOXER_MAIN_SOURCE_FILES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}" "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp")

foreach (main_file ${FLOXER_MAIN_SOURCE_FILES})
    get_filename_component (target_name ${main_file} NAME_WE)

    if (target_name STREQUAL "${PROJECT_NAME}")
        add_executable (${target_name} ${main_file})
    else ()
        add_executable (${target_name} EXCLUDE_FROM_ALL ${main_file})
    endif()

    target_link_libraries (${target_name} PUBLIC "${PROJECT_NAME}_lib")

    unset (main_file)
    unset (target_name)
endforeach ()
