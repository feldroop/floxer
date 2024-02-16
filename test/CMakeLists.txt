CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    GIT_TAG release-1.12.1
    VERSION 1.12.1
    OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
)

add_executable(search_test search_test.cpp)
target_link_libraries(
    search_test
    GTest::gtest_main
    pex_lib
    search_lib
)

add_executable(pex_tree_test pex_tree_test.cpp)
target_link_libraries(
    pex_tree_test
    GTest::gtest_main
    pex_lib
)

add_test(NAME search COMMAND search_test)
add_test(NAME pex_tree COMMAND pex_tree_test)
