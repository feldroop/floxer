cmake_minimum_required (VERSION 3.16)

project (floxer VERSION 0.2.0)

# Force language standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Use ccache.
include (cmake/ccache.cmake)
require_ccache ()

# for the CPM package manager
include(cmake/CPM.cmake)

# dependencies
CPMaddPackage(
    NAME cereal
    GITHUB_REPOSITORY USCiLab/cereal
    GIT_TAG v1.3.2
    OPTIONS
        "SKIP_PERFORMANCE_COMPARISON ON"
        "BUILD_SANDBOX OFF"
        "BUILD_DOC OFF"
        "WITH_WERROR OFF"
    EXCLUDE_FROM_ALL TRUE
    SYSTEM YES
)
find_package(sharg REQUIRED PATHS lib/sharg/build_system)
find_package(ivio REQUIRED PATHS lib/ivio)
find_package(ivsigma REQUIRED PATHS lib/ivsigma)
find_package(fmindex_collection REQUIRED PATHS lib/fmindex)

# floxer source code
add_subdirectory(src)

# testing
list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
