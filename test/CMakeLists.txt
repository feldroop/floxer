CPMGetPackage(googletest)

list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
enable_testing ()

include(GoogleTest)

add_custom_target (check COMMAND ${CMAKE_CTEST_COMMAND} ${CMAKE_CTEST_ARGUMENTS})
file(GLOB FLOXER_TEST_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)

# loop through all files in this folder and declare a test for each of them
foreach (test_file ${FLOXER_TEST_FILES})
    get_filename_component (test_name ${test_file} NAME_WE)

    add_executable (${test_name} ${test_file})
    target_link_libraries (${test_name} PUBLIC "${PROJECT_NAME}_lib" GTest::gtest_main)
    target_compile_options (${test_name} PUBLIC "-pedantic" "-Wall" "-Wextra")

    add_dependencies (${test_name} ${PROJECT_NAME})
    add_dependencies (check ${test_name})

    gtest_discover_tests (${test_name})

    unset (test_name)
    unset (test_file)
endforeach ()

message (STATUS "You can run `make check` to build and run tests.")
