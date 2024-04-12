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

# Use ccache for faster cached compiles
include (cmake/ccache.cmake)
require_ccache ()

# for the CPM package manager
include(cmake/CPM.cmake)

# ----- dependencies -----

# for FM-index serialization
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

# for logging and output (includes libfmt)
CPMaddPackage("gh:gabime/spdlog#v1.13.0")

# for CLI parsing
CPMaddPackage(
    NAME sharg
    GITHUB_REPOSITORY seqan/sharg-parser
    GIT_TAG 39f65a4890f8c5108af2b5c7974893ff6ed87e50
    OPTIONS
        "SHARG_NO_TDL ON"
    EXCLUDE_FROM_ALL TRUE
    SYSTEM YES
)

# for sequence input
CPMaddPackage("gh:iv-project/IVio#8fcb901bc08dd48d050bcfbcf73d5fe40db6c90a")

# for sequence alphabet and transformation utilities
CPMaddPackage("gh:iv-project/IVSigma#4b76036ed331b91b7a08fc7249ac51b4def71c74")

# for the FM-index
# set(LIBSAIS_USE_OPENMP ON)
CPMaddPackage("gh:SGSSGene/fmindex-collection#b0e311f122f6ae9d15d470ca88ac21f028dc010a")

# ----- floxer source code -----

add_subdirectory(src)

# testing
list (APPEND CMAKE_CTEST_ARGUMENTS "--output-on-failure")
enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
