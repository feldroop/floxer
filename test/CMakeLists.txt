CPMAddPackage(
    NAME googletest
    GITHUB_REPOSITORY google/googletest
    GIT_TAG release-1.12.1
    VERSION 1.12.1
    OPTIONS "INSTALL_GTEST OFF" "gtest_force_shared_crt"
)

add_executable(verification_test verification_test.cpp)
target_link_libraries(
    verification_test
    GTest::gtest_main
    verification_lib
)

add_test(NAME verification COMMAND verification_test)
# TODO add more tests once the core of the program is more or less stable
